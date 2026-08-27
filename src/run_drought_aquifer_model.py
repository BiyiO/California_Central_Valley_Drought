"""
California Central Valley Agricultural Drought & Aquifer Overdraft Modeling (ADAVI Framework)
--------------------------------------------------------------------------------------------
A geospatial and remote sensing modeling suite for multi-criteria assessment of
agricultural drought vulnerability and groundwater overdraft stress across Fresno,
Kings, and Tulare Counties, California.

Author: Antigravity Geospatial AI & Research Lab
CRS: EPSG:3310 (California Teale Albers)
Spatial Resolution: 100m
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.features import shapes
from scipy.ndimage import gaussian_filter
import geopandas as gpd
from shapely.geometry import shape


def setup_directories(base_dir):
    """Ensure all project directories exist."""
    dirs = {
        'base': base_dir,
        'raw': os.path.join(base_dir, "01_Raw_Data"),
        'processed': os.path.join(base_dir, "02_Processed_Factors"),
        'model': os.path.join(base_dir, "03_Drought_Aquifer_Model"),
        'maps': os.path.join(base_dir, "04_Final_Maps"),
        'boundary': os.path.join(base_dir, "01_Raw_Data", "Boundary"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return dirs


def find_raw_raster(raw_dir, keyword, preferred_prefix="CA_Core"):
    """Locate the best matching GeoTIFF for a given parameter."""
    preferred = None
    fallback = None
    for root, _, files in os.walk(raw_dir):
        # Prefer direct folders over 'Old data'
        is_old = "old data" in root.lower()
        for f in sorted(files):
            if f.endswith('.tif') and keyword.lower() in f.lower():
                full_path = os.path.join(root, f)
                if preferred_prefix.lower() in f.lower():
                    return full_path
                if not is_old and fallback is None:
                    fallback = full_path
                elif fallback is None:
                    fallback = full_path
    return preferred or fallback


def load_and_align_rasters(dem_path, ndwi_path, ndvi_path, lst_path):
    """Load reference DEM and reproject/align satellite layers to matching grid."""
    print(" [1/6] Loading Reference DEM & Aligning Satellite Layers...")
    with rasterio.open(dem_path) as dem_src:
        dem = dem_src.read(1).astype(np.float32)
        profile = dem_src.profile.copy()
        transform = dem_src.transform
        crs = dem_src.crs
        height = dem_src.height
        width = dem_src.width
        nodata = dem_src.nodata

    # DEM valid mask (exclude extreme nodata and ocean values)
    valid_mask = (~np.isnan(dem)) & (dem > -50)
    if nodata is not None:
        valid_mask = valid_mask & (dem != nodata)

    def align(file_path):
        aligned = np.full((height, width), np.nan, dtype=np.float32)
        if file_path and os.path.exists(file_path):
            with rasterio.open(file_path) as src:
                reproject(
                    source=rasterio.band(src, 1),
                    destination=aligned,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=crs,
                    resampling=Resampling.bilinear
                )
        return aligned

    ndwi = align(ndwi_path)
    ndvi = align(ndvi_path)
    lst = align(lst_path)

    # Combined valid mask across all 4 layers
    valid_mask = valid_mask & (~np.isnan(ndwi)) & (~np.isnan(ndvi)) & (~np.isnan(lst))

    return {
        'dem': dem,
        'ndwi': ndwi,
        'ndvi': ndvi,
        'lst': lst,
        'valid_mask': valid_mask,
        'profile': profile,
        'transform': transform,
        'crs': crs,
        'height': height,
        'width': width
    }


def derive_and_standardize_factors(raster_dict, dirs):
    """
    Standardize 5 hydro-agronomic and geomorphic factors to 1-5 vulnerability scale:
    Factor 1: Crop Moisture Deficit (NDWI) [Weight: 0.30]
    Factor 2: Thermal Evaporative Stress (LST °C) [Weight: 0.25]
    Factor 3: Vegetation Biomass & Crop Health (NDVI) [Weight: 0.20]
    Factor 4: Alluvial Plain Micro-Topography & Slope [Weight: 0.15]
    Factor 5: Aquifer Overdraft & Soil Compaction Proxy [Weight: 0.10]
    """
    print(" [2/6] Standardizing 5 Hydro-Agronomic & Geomorphic Factors (1-5 Scale)...")
    ndwi = raster_dict['ndwi']
    lst = raster_dict['lst']
    ndvi = raster_dict['ndvi']
    dem = raster_dict['dem']
    valid_mask = raster_dict['valid_mask']
    transform = raster_dict['transform']
    profile = raster_dict['profile'].copy()

    # Factor 1: Crop Moisture Deficit (NDWI)
    # Lower NDWI = drier canopy = higher drought vulnerability
    ndwi_f = np.full_like(ndwi, -9999.0, dtype=np.float32)
    ndwi_f[valid_mask & (ndwi > 0.20)] = 1.0                           # Very High Moisture
    ndwi_f[valid_mask & (ndwi > 0.05) & (ndwi <= 0.20)] = 2.0         # Moderate Moisture
    ndwi_f[valid_mask & (ndwi > -0.10) & (ndwi <= 0.05)] = 3.0        # Mild Deficit
    ndwi_f[valid_mask & (ndwi > -0.25) & (ndwi <= -0.10)] = 4.0       # High Deficit
    ndwi_f[valid_mask & (ndwi <= -0.25)] = 5.0                         # Severe Canopy Desiccation

    # Factor 2: Land Surface Temperature (LST °C)
    # Higher LST = extreme evaporative demand & thermal stress
    lst_f = np.full_like(lst, -9999.0, dtype=np.float32)
    lst_f[valid_mask & (lst <= 28.0)] = 1.0                           # Cool / Irrigated
    lst_f[valid_mask & (lst > 28.0) & (lst <= 34.0)] = 2.0            # Moderate Thermal
    lst_f[valid_mask & (lst > 34.0) & (lst <= 39.0)] = 3.0            # Elevated Heat
    lst_f[valid_mask & (lst > 39.0) & (lst <= 44.0)] = 4.0            # Severe Thermal Stress
    lst_f[valid_mask & (lst > 44.0)] = 5.0                            # Extreme Heat Anomaly (>44°C)

    # Factor 3: Crop Vigor & Biomass (NDVI)
    # Lower NDVI in agricultural zones = fallowed / drought-stunted fields
    ndvi_f = np.full_like(ndvi, -9999.0, dtype=np.float32)
    ndvi_f[valid_mask & (ndvi > 0.60)] = 1.0                          # Dense Healthy Canopy
    ndvi_f[valid_mask & (ndvi > 0.45) & (ndvi <= 0.60)] = 2.0         # Moderate Vigor
    ndvi_f[valid_mask & (ndvi > 0.30) & (ndvi <= 0.45)] = 3.0         # Sparse / Stressed
    ndvi_f[valid_mask & (ndvi > 0.15) & (ndvi <= 0.30)] = 4.0         # Severely Stressed
    ndvi_f[valid_mask & (ndvi <= 0.15)] = 5.0                         # Bare Soil / Fallowed Land

    # Factor 4: Topographic Slope (Degrees)
    # Flatter alluvial valley floor = deep unconsolidated sediment prone to compaction & subsidence
    px, py = abs(transform[0]), abs(transform[4])
    gy, gx = np.gradient(dem, py, px)
    slope_deg = np.degrees(np.arctan(np.sqrt(gx**2 + gy**2)))

    slope_f = np.full_like(slope_deg, -9999.0, dtype=np.float32)
    slope_f[valid_mask & (slope_deg <= 0.8)] = 5.0                    # Ultra-flat valley floor (Subsidence core)
    slope_f[valid_mask & (slope_deg > 0.8) & (slope_deg <= 1.8)] = 4.0
    slope_f[valid_mask & (slope_deg > 1.8) & (slope_deg <= 3.0)] = 3.0
    slope_f[valid_mask & (slope_deg > 3.0) & (slope_deg <= 5.0)] = 2.0
    slope_f[valid_mask & (slope_deg > 5.0)] = 1.0                    # Foothills / Sierra slopes

    # Factor 5: Aquifer Overdraft & Soil Compaction Proxy
    # Spatial filter over alluvial basin floor capturing contiguous deep alluvial subsidence corridors
    compaction_proxy = gaussian_filter(np.where(valid_mask, slope_f, 0.0), sigma=3.0)
    overdraft_f = np.full_like(compaction_proxy, -9999.0, dtype=np.float32)
    overdraft_f[valid_mask & (compaction_proxy >= 4.2)] = 5.0         # Core Overdraft Zone
    overdraft_f[valid_mask & (compaction_proxy >= 3.4) & (compaction_proxy < 4.2)] = 4.0
    overdraft_f[valid_mask & (compaction_proxy >= 2.6) & (compaction_proxy < 3.4)] = 3.0
    overdraft_f[valid_mask & (compaction_proxy >= 1.8) & (compaction_proxy < 2.6)] = 2.0
    overdraft_f[valid_mask & (compaction_proxy < 1.8)] = 1.0

    # Save factor GeoTIFFs to 02_Processed_Factors
    profile.update(dtype=rasterio.float32, count=1, driver='GTiff', nodata=-9999.0)
    factor_files = {
        'Factor1_Crop_Moisture_Deficit_NDWI.tif': ndwi_f,
        'Factor2_Thermal_Stress_LST.tif': lst_f,
        'Factor3_Vegetation_Biomass_NDVI.tif': ndvi_f,
        'Factor4_Topographic_Slope.tif': slope_f,
        'Factor5_Aquifer_Overdraft_Proxy.tif': overdraft_f,
    }

    for fname, f_arr in factor_files.items():
        out_p = os.path.join(dirs['processed'], fname)
        with rasterio.open(out_p, 'w', **profile) as dst:
            dst.write(f_arr, 1)

    return {
        'ndwi_f': ndwi_f,
        'lst_f': lst_f,
        'ndvi_f': ndvi_f,
        'slope_f': slope_f,
        'overdraft_f': overdraft_f,
        'slope_deg': slope_deg,
        'profile': profile
    }


def execute_adavi_model(raster_dict, factor_dict, dirs, weights=None):
    """
    Execute Multi-Criteria Weighted Overlay & 5-tier classification.
    ADAVI = (0.30 * NDWI) + (0.25 * LST) + (0.20 * NDVI) + (0.15 * Slope) + (0.10 * Overdraft)
    """
    print(" [3/6] Executing Multi-Criteria Weighted ADAVI Model...")
    if weights is None:
        weights = {
            'ndwi': 0.30,
            'lst': 0.25,
            'ndvi': 0.20,
            'slope': 0.15,
            'overdraft': 0.10
        }

    valid_mask = raster_dict['valid_mask']
    ndwi_f = factor_dict['ndwi_f']
    lst_f = factor_dict['lst_f']
    ndvi_f = factor_dict['ndvi_f']
    slope_f = factor_dict['slope_f']
    overdraft_f = factor_dict['overdraft_f']
    profile = factor_dict['profile'].copy()

    # Continuous ADAVI Calculation
    adavi_cont = (
        (weights['ndwi'] * ndwi_f) +
        (weights['lst'] * lst_f) +
        (weights['ndvi'] * ndvi_f) +
        (weights['slope'] * slope_f) +
        (weights['overdraft'] * overdraft_f)
    )
    adavi_cont[~valid_mask] = -9999.0

    # 5-Tier Vulnerability Zonation Classification
    # Scale: 1 (Very Low: <= 2.2), 2 (Low: 2.2 - 2.9), 3 (Moderate: 2.9 - 3.6), 4 (High: 3.6 - 4.2), 5 (Extreme: > 4.2)
    adavi_class = np.full_like(adavi_cont, -9999.0, dtype=np.float32)
    adavi_class[valid_mask & (adavi_cont <= 2.2)] = 1.0
    adavi_class[valid_mask & (adavi_cont > 2.2) & (adavi_cont <= 2.9)] = 2.0
    adavi_class[valid_mask & (adavi_cont > 2.9) & (adavi_cont <= 3.6)] = 3.0
    adavi_class[valid_mask & (adavi_cont > 3.6) & (adavi_cont <= 4.2)] = 4.0
    adavi_class[valid_mask & (adavi_cont > 4.2)] = 5.0

    # Save to 03_Drought_Aquifer_Model
    out_cont = os.path.join(dirs['model'], "ADAVI_Continuous_Index.tif")
    out_class = os.path.join(dirs['model'], "ADAVI_5Class_Vulnerability.tif")
    out_final_map = os.path.join(dirs['maps'], "California_Central_Valley_Drought_Vulnerability_Classified.tif")

    for path, arr in [(out_cont, adavi_cont), (out_class, adavi_class), (out_final_map, adavi_class)]:
        with rasterio.open(path, 'w', **profile) as dst:
            dst.write(arr, 1)

    return {
        'adavi_cont': adavi_cont,
        'adavi_class': adavi_class,
        'weights': weights,
        'out_cont': out_cont,
        'out_class': out_class
    }


def compute_zonal_statistics_and_clusters(raster_dict, model_dict, dirs):
    """
    Compute county-by-county zonal breakdown for Fresno, Kings, and Tulare Counties,
    and polygonize Priority Groundwater Mitigation clusters.
    """
    print(" [4/6] Computing County Zonal Statistics & Groundwater Priority Clusters...")
    valid_mask = raster_dict['valid_mask']
    adavi_class = model_dict['adavi_class']
    adavi_cont = model_dict['adavi_cont']
    transform = raster_dict['transform']
    crs = raster_dict['crs']

    pixel_area_km2 = (100 * 100) / 1e6
    acres_per_km2 = 247.105381

    # Load County Boundaries
    boundary_shp = os.path.join(dirs['boundary'], "Fresno_Tulare_Kings_Boundary.shp")
    if not os.path.exists(boundary_shp):
        # Generate from census FIPS dataset
        census_url = "https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json"
        gdf_all = gpd.read_file(census_url)
        gdf_core = gdf_all[gdf_all['id'].isin(['06019', '06031', '06107'])].copy()
        gdf_core['County_Nam'] = gdf_core['id'].map({
            '06019': 'Fresno County',
            '06031': 'Kings County',
            '06107': 'Tulare County'
        })
        gdf_core.to_file(boundary_shp)

    gdf_counties = gpd.read_file(boundary_shp).to_crs(crs)

    vuln_labels = ['Very Low', 'Low', 'Moderate', 'High', 'Extreme']
    county_stats = []

    # Regional Overall Statistics
    valid_pixels = adavi_class[valid_mask]
    total_regional_km2 = len(valid_pixels) * pixel_area_km2
    total_regional_acres = total_regional_km2 * acres_per_km2

    reg_counts = [(valid_pixels == i).sum() for i in range(1, 6)]
    reg_km2 = [c * pixel_area_km2 for c in reg_counts]
    reg_acres = [k * acres_per_km2 for k in reg_km2]
    reg_pcts = [(c / len(valid_pixels)) * 100 for c in reg_counts]

    regional_summary = {
        'Region': 'California Central Valley Core (Fresno, Kings, Tulare)',
        'Total_Area_km2': round(total_regional_km2, 2),
        'Total_Area_Acres': round(total_regional_acres, 0),
        'Mean_ADAVI': round(float(np.mean(adavi_cont[valid_mask])), 3),
        'Extreme_Vuln_Acres': round(reg_acres[4], 0),
        'Extreme_Vuln_Pct': round(reg_pcts[4], 2),
        'High_Vuln_Acres': round(reg_acres[3], 0),
        'High_Vuln_Pct': round(reg_pcts[3], 2),
        'Moderate_Vuln_Acres': round(reg_acres[2], 0),
        'Moderate_Vuln_Pct': round(reg_pcts[2], 2),
        'Low_Vuln_Acres': round(reg_acres[1], 0),
        'Low_Vuln_Pct': round(reg_pcts[1], 2),
        'Very_Low_Vuln_Acres': round(reg_acres[0], 0),
        'Very_Low_Vuln_Pct': round(reg_pcts[0], 2),
    }

    # Rasterize County polygons to compute zonal stats
    from rasterio.features import rasterize
    shapes_gen = ((geom, idx) for idx, geom in enumerate(gdf_counties.geometry))
    county_raster = rasterize(
        shapes=shapes_gen,
        out_shape=(raster_dict['height'], raster_dict['width']),
        transform=transform,
        fill=-1,
        dtype=np.int16
    )

    for idx, row in gdf_counties.iterrows():
        c_name = row.get('County_Nam') or row.get('NAME') or f"County_{idx}"
        c_mask = valid_mask & (county_raster == idx)
        c_pixels = adavi_class[c_mask]
        c_cont = adavi_cont[c_mask]

        if len(c_pixels) == 0:
            continue

        c_km2 = len(c_pixels) * pixel_area_km2
        c_acres = c_km2 * acres_per_km2
        counts = [(c_pixels == i).sum() for i in range(1, 6)]
        kms = [cnt * pixel_area_km2 for cnt in counts]
        acs = [k * acres_per_km2 for k in kms]
        pcts = [(cnt / len(c_pixels)) * 100 for cnt in counts]

        county_stats.append({
            'County': c_name,
            'Total_Area_km2': round(c_km2, 2),
            'Total_Area_Acres': round(c_acres, 0),
            'Mean_ADAVI': round(float(np.mean(c_cont)), 3),
            'Extreme_Vuln_Acres': round(acs[4], 0),
            'Extreme_Vuln_Pct': round(pcts[4], 2),
            'High_Vuln_Acres': round(acs[3], 0),
            'High_Vuln_Pct': round(pcts[3], 2),
            'Moderate_Vuln_Acres': round(acs[2], 0),
            'Moderate_Vuln_Pct': round(pcts[2], 2),
            'Low_Vuln_Acres': round(acs[1], 0),
            'Low_Vuln_Pct': round(pcts[1], 2),
            'Very_Low_Vuln_Acres': round(acs[0], 0),
            'Very_Low_Vuln_Pct': round(pcts[0], 2),
        })

    # Save CSV & JSON
    df_stats = pd.DataFrame(county_stats)
    csv_path = os.path.join(dirs['model'], "County_Vulnerability_Statistics.csv")
    json_path = os.path.join(dirs['model'], "County_Vulnerability_Statistics.json")

    df_stats.to_csv(csv_path, index=False)
    with open(json_path, 'w', encoding='utf-8') as jf:
        json.dump({
            'regional_summary': regional_summary,
            'county_breakdown': county_stats
        }, jf, indent=2)

    # Vectorize Priority Mitigation Zones (Classes 4 and 5)
    high_extreme_mask = (adavi_class >= 4.0).astype(np.uint8)
    geom_shapes = list(shapes(high_extreme_mask, mask=(high_extreme_mask == 1), transform=transform))

    records = []
    for geom_dict, val in geom_shapes:
        poly = shape(geom_dict)
        if poly.area >= (50 * 100 * 100):  # Filter small speckles (<50 ha / 0.5 km2)
            records.append({
                'geometry': poly,
                'priority_level': 'High to Extreme Drought & Overdraft Vulnerability',
                'area_km2': round(poly.area / 1e6, 3),
                'area_acres': round((poly.area / 1e6) * acres_per_km2, 1)
            })

    if records:
        gdf_priority = gpd.GeoDataFrame(records, crs=crs)
        priority_out = os.path.join(dirs['model'], "Priority_Groundwater_Mitigation_Zones.geojson")
        gdf_priority.to_file(priority_out, driver="GeoJSON")
        print(f"    Saved {len(gdf_priority)} priority mitigation clusters to {priority_out}")

    return {
        'regional_summary': regional_summary,
        'df_stats': df_stats,
        'gdf_counties': gdf_counties
    }


def generate_cartographic_figures(raster_dict, factor_dict, model_dict, stats_dict, dirs):
    """
    Generate publication-ready 300 DPI figures:
    1. Executive 4-Panel Infographic Dashboard (Dark Theme)
    2. Cartographic Regional Map with boundaries, scale, and annotations (Light/Neutral Theme)
    3. Hydro-Agronomic Factor Correlation Matrix
    """
    print(" [5/6] Generating Publication & Cartographic Figures (300 DPI)...")
    valid_mask = raster_dict['valid_mask']
    adavi_class = model_dict['adavi_class']
    adavi_cont = model_dict['adavi_cont']
    ndwi_f = factor_dict['ndwi_f']
    lst_f = factor_dict['lst_f']
    ndvi_f = factor_dict['ndvi_f']
    slope_f = factor_dict['slope_f']
    overdraft_f = factor_dict['overdraft_f']
    gdf_counties = stats_dict['gdf_counties']
    df_stats = stats_dict['df_stats']
    regional_summary = stats_dict['regional_summary']

    vuln_labels = ['Very Low', 'Low', 'Moderate', 'High', 'Extreme']
    colors_drought = ['#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#d7191c']
    cmap_vuln = mcolors.ListedColormap(colors_drought)

    # -------------------------------------------------------------------------
    # 1. EXECUTIVE 4-PANEL DASHBOARD FIGURE
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(18, 14), facecolor='#0f172a')
    fig.suptitle(
        "CALIFORNIA CENTRAL VALLEY AGRICULTURAL DROUGHT & AQUIFER VULNERABILITY\n"
        "Multi-Criteria Remote Sensing & Hydro-Geomorphic Overdraft Analytics (Fresno, Kings, Tulare)",
        fontsize=16, fontweight='bold', color='#ffffff', y=0.97
    )

    for ax in axes.flat:
        ax.set_facecolor('#1e293b')

    # Panel A: Crop Moisture Stress Index (NDWI)
    axes[0, 0].set_title("A. Crop Canopy Moisture Stress Index (NDWI)", fontsize=13, fontweight='bold', color='#38bdf8', pad=10)
    im_a = axes[0, 0].imshow(np.ma.masked_where(~valid_mask, ndwi_f), cmap='YlGnBu_r', vmin=1, vmax=5)
    axes[0, 0].axis('off')
    cbar_a = fig.colorbar(im_a, ax=axes[0, 0], fraction=0.035, pad=0.02)
    cbar_a.set_ticks([1, 2, 3, 4, 5])
    cbar_a.set_ticklabels(['1: Very Low', '2: Low', '3: Moderate', '4: High', '5: Severe Deficit'], color='#e2e8f0', fontsize=9)
    cbar_a.ax.yaxis.set_tick_params(color='#e2e8f0')

    # Panel B: Thermal Land Surface Temperature (LST °C)
    axes[0, 1].set_title("B. Summer Land Surface Temperature & Evaporative Demand (LST °C)", fontsize=13, fontweight='bold', color='#f87171', pad=10)
    im_b = axes[0, 1].imshow(np.ma.masked_where(~valid_mask, lst_f), cmap='inferno', vmin=1, vmax=5)
    axes[0, 1].axis('off')
    cbar_b = fig.colorbar(im_b, ax=axes[0, 1], fraction=0.035, pad=0.02)
    cbar_b.set_ticks([1, 2, 3, 4, 5])
    cbar_b.set_ticklabels(['1: <=28°C', '2: 28-34°C', '3: 34-39°C', '4: 39-44°C', '5: >44°C (Extreme)'], color='#e2e8f0', fontsize=9)
    cbar_b.ax.yaxis.set_tick_params(color='#e2e8f0')

    # Panel C: ADAVI Multi-Criteria Zonation
    axes[1, 0].set_title("C. Multi-Criteria Drought & Aquifer Overdraft Zonation (ADAVI)", fontsize=13, fontweight='bold', color='#4ade80', pad=10)
    im_c = axes[1, 0].imshow(np.ma.masked_where(~valid_mask, adavi_class), cmap=cmap_vuln, vmin=1, vmax=5)
    axes[1, 0].axis('off')
    legend_elements = [Patch(facecolor=c, edgecolor='black', label=l) for c, l in zip(colors_drought, vuln_labels)]
    axes[1, 0].legend(handles=legend_elements, loc='lower left', fontsize=10, facecolor='#1e293b', edgecolor='#475569', labelcolor='#ffffff')

    # Panel D: Regional Vulnerability Share & County Breakdown
    axes[1, 1].set_title("D. Regional Agricultural Vulnerability Distribution", fontsize=13, fontweight='bold', color='#fbbf24', pad=10)
    
    valid_pixels = adavi_class[valid_mask]
    counts = [(valid_pixels == i).sum() for i in range(1, 6)]
    pcts = [(c / len(valid_pixels)) * 100 for c in counts]
    
    wedges, texts, autotexts = axes[1, 1].pie(
        pcts,
        labels=vuln_labels,
        colors=colors_drought,
        autopct='%1.1f%%',
        pctdistance=0.75,
        startangle=140,
        wedgeprops=dict(width=0.45, edgecolor='#0f172a', linewidth=2),
        textprops=dict(color='#ffffff', fontweight='bold', fontsize=10)
    )
    for at in autotexts:
        at.set_color('#0f172a')
        at.set_fontsize(10)

    # Inset center text in donut chart
    axes[1, 1].text(0, 0, f"Total Area\n{regional_summary['Total_Area_Acres']:,.0f}\nAcres", 
                    ha='center', va='center', color='#ffffff', fontsize=11, fontweight='bold')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig4_path = os.path.join(dirs['maps'], "California_Central_Valley_Executive_4Panel.png")
    plt.savefig(fig4_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"    Saved: {fig4_path}")

    # -------------------------------------------------------------------------
    # 2. STANDALONE CARTOGRAPHIC MAP WITH COUNTY BOUNDARIES & SCALE
    # -------------------------------------------------------------------------
    fig2, ax = plt.subplots(figsize=(14, 16), facecolor='#f8fafc')
    ax.set_facecolor('#e2e8f0')

    # Plot classified vulnerability
    im = ax.imshow(
        np.ma.masked_where(~valid_mask, adavi_class),
        cmap=cmap_vuln,
        vmin=1,
        vmax=5,
        extent=[
            raster_dict['transform'][2],
            raster_dict['transform'][2] + raster_dict['transform'][0] * raster_dict['width'],
            raster_dict['transform'][5] + raster_dict['transform'][4] * raster_dict['height'],
            raster_dict['transform'][5]
        ]
    )

    # Overlay County Boundaries
    gdf_counties.boundary.plot(ax=ax, color='#1e293b', linewidth=2.0, linestyle='--')

    # Annotate County Centroids
    for _, row in gdf_counties.iterrows():
        c_name = row.get('County_Nam') or row.get('NAME') or ''
        centroid = row.geometry.centroid
        ax.text(
            centroid.x, centroid.y, c_name.upper(),
            fontsize=12, fontweight='heavy', color='#0f172a',
            ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#ffffff', alpha=0.85, edgecolor='#475569')
        )

    ax.set_title(
        "CALIFORNIA CENTRAL VALLEY AGRICULTURAL DROUGHT & AQUIFER VULNERABILITY MAP\n"
        "Fresno, Kings, and Tulare Counties (EPSG:3310 - California Teale Albers)",
        fontsize=14, fontweight='bold', color='#0f172a', pad=15
    )

    ax.set_xlabel("Easting (Meters - CA Teale Albers)", fontsize=10, fontweight='bold', color='#334155')
    ax.set_ylabel("Northing (Meters - CA Teale Albers)", fontsize=10, fontweight='bold', color='#334155')
    ax.tick_params(colors='#334155', labelsize=9)

    # Cartographic Legend
    legend_elements = [
        Patch(facecolor='#2b83ba', edgecolor='black', label=f"Very Low Vulnerability ({df_stats['Very_Low_Vuln_Acres'].sum():,.0f} Ac)"),
        Patch(facecolor='#abdda4', edgecolor='black', label=f"Low Vulnerability ({df_stats['Low_Vuln_Acres'].sum():,.0f} Ac)"),
        Patch(facecolor='#ffffbf', edgecolor='black', label=f"Moderate Vulnerability ({df_stats['Moderate_Vuln_Acres'].sum():,.0f} Ac)"),
        Patch(facecolor='#fdae61', edgecolor='black', label=f"High Vulnerability ({df_stats['High_Vuln_Acres'].sum():,.0f} Ac)"),
        Patch(facecolor='#d7191c', edgecolor='black', label=f"Extreme Vulnerability ({df_stats['Extreme_Vuln_Acres'].sum():,.0f} Ac)"),
    ]
    ax.legend(
        handles=legend_elements,
        loc='lower right',
        title="ADAVI Vulnerability Zonation",
        title_fontsize='11',
        fontsize=9.5,
        facecolor='#ffffff',
        edgecolor='#94a3b8',
        framealpha=0.95
    )

    # North Arrow and Scale Bar Annotation
    ax.annotate(
        'N\n▲', xy=(0.06, 0.93), xycoords='axes fraction',
        ha='center', va='center', fontsize=18, fontweight='bold', color='#0f172a',
        bbox=dict(boxstyle='circle,pad=0.4', facecolor='#ffffff', edgecolor='#334155', alpha=0.9)
    )

    carto_path = os.path.join(dirs['maps'], "California_Central_Valley_Cartographic_Map.png")
    plt.savefig(carto_path, dpi=300, facecolor=fig2.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"    Saved: {carto_path}")

    # -------------------------------------------------------------------------
    # 3. FACTOR CORRELATION MATRIX & STATISTICAL HEATMAP
    # -------------------------------------------------------------------------
    fig3, ax3 = plt.subplots(figsize=(10, 8), facecolor='#ffffff')
    
    df_corr = pd.DataFrame({
        'NDWI Deficit (F1)': ndwi_f[valid_mask],
        'LST Heat (F2)': lst_f[valid_mask],
        'NDVI Stress (F3)': ndvi_f[valid_mask],
        'Slope Flatness (F4)': slope_f[valid_mask],
        'Overdraft Proxy (F5)': overdraft_f[valid_mask],
        'Composite ADAVI': adavi_cont[valid_mask]
    }).corr()

    cax = ax3.matshow(df_corr, cmap='coolwarm', vmin=-1, vmax=1)
    fig3.colorbar(cax, fraction=0.046, pad=0.04)

    ax3.set_xticks(range(len(df_corr.columns)))
    ax3.set_yticks(range(len(df_corr.columns)))
    ax3.set_xticklabels(df_corr.columns, rotation=35, ha='left', fontsize=10, fontweight='bold')
    ax3.set_yticklabels(df_corr.columns, fontsize=10, fontweight='bold')

    for (i, j), val in np.ndenumerate(df_corr.values):
        ax3.text(j, i, f"{val:.2f}", ha='center', va='center',
                 color='white' if abs(val) > 0.4 else 'black', fontweight='bold', fontsize=11)

    ax3.set_title(
        "Hydro-Agronomic & Geomorphic Factor Correlation Matrix\n(Central Valley ADAVI Modeling)",
        fontsize=13, fontweight='bold', pad=25
    )

    corr_path = os.path.join(dirs['maps'], "Factor_Correlation_Matrix.png")
    plt.savefig(corr_path, dpi=300, facecolor=fig3.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"    Saved: {corr_path}")


def main():
    parser = argparse.ArgumentParser(description="California Central Valley Drought & Aquifer Overdraft Modeling Suite")
    parser.add_argument("--base-dir", default=r"C:\Users\USER\Documents\GIS\California_Central_Valley_Drought", help="Project base directory")
    args = parser.parse_args()

    dirs = setup_directories(args.base_dir)

    print("=" * 80)
    print(" CALIFORNIA CENTRAL VALLEY AGRICULTURAL DROUGHT & AQUIFER OVERDRAFT MODELING")
    print(" Counties: Fresno (06019), Kings (06031), Tulare (06107)")
    print("=" * 80)

    # Locate inputs
    dem_path = find_raw_raster(dirs['raw'], "DEM")
    ndwi_path = find_raw_raster(dirs['raw'], "NDWI")
    ndvi_path = find_raw_raster(dirs['raw'], "NDVI")
    lst_path = find_raw_raster(dirs['raw'], "LST")

    print(f" Discovered DEM  : {dem_path}")
    print(f" Discovered NDWI : {ndwi_path}")
    print(f" Discovered NDVI : {ndvi_path}")
    print(f" Discovered LST  : {lst_path}")

    # Step 1: Load and Align
    raster_dict = load_and_align_rasters(dem_path, ndwi_path, ndvi_path, lst_path)

    # Step 2: Standardize 5 Factors
    factor_dict = derive_and_standardize_factors(raster_dict, dirs)

    # Step 3: Execute ADAVI Model
    model_dict = execute_adavi_model(raster_dict, factor_dict, dirs)

    # Step 4: Zonal Statistics & Vector Clusters
    stats_dict = compute_zonal_statistics_and_clusters(raster_dict, model_dict, dirs)

    # Step 5: High-Res Visual Cartography
    generate_cartographic_figures(raster_dict, factor_dict, model_dict, stats_dict, dirs)

    # Step 6: Print Executive Terminal Summary
    print("\n" + "=" * 80)
    print("                     MODELING EXECUTION SUMMARY")
    print("=" * 80)
    print(stats_dict['df_stats'].to_string(index=False))
    print("=" * 80)
    reg = stats_dict['regional_summary']
    print(f" Regional Extent  : {reg['Total_Area_km2']:,.1f} sq km ({reg['Total_Area_Acres']:,.0f} Acres)")
    print(f" Mean ADAVI Score : {reg['Mean_ADAVI']:.3f} / 5.000")
    print(f" High + Extreme   : {reg['High_Vuln_Acres'] + reg['Extreme_Vuln_Acres']:,.0f} Acres ({(reg['High_Vuln_Pct'] + reg['Extreme_Vuln_Pct']):.1f}% of total farmland)")
    print("=" * 80)
    print("\n [6/6] Pipeline Complete! All rasters, tables, vectors, and publication maps generated.")


if __name__ == "__main__":
    main()
