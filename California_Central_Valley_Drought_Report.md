# California Central Valley Agricultural Drought & Aquifer Overdraft Modeling
## Geospatial Decision Support & Vulnerability Zonation Report (ADAVI Framework)

**Focus Area:** Fresno, Kings, and Tulare Counties, California  
**Geographic Domain:** San Joaquin Valley & Tulare Lake Hydrologic Basins  
**Spatial Resolution:** 100-Meter Grid  
**Projected Coordinate System:** EPSG:3310 (California Teale Albers, NAD83)  
**Analysis Epoch:** Peak Summer Agricultural Irrigation Season (Landsat 8/9 OLI/TIRS & USGS 3DEP Analytics)  

---

## 1. Executive Summary

California's Central Valley—specifically the three-county core comprising **Fresno, Kings, and Tulare Counties**—produces a substantial portion of the nation's fruit, nut, and dairy commodities. However, intensifying hydro-climatic volatility, chronic groundwater overdraft, deep alluvial soil compaction, and surface water delivery curtailments threaten long-term agricultural sustainability.

To provide spatial clarity for **Groundwater Sustainability Agencies (GSAs)**, regional water managers, and agricultural planners under the **Sustainable Groundwater Management Act (SGMA)**, this study deploys the **Agricultural Drought & Aquifer Vulnerability Index (ADAVI)**. The ADAVI model integrates five remote sensing and hydro-geomorphic factors across **7,848,462 acres (31,761.6 km²)** of land.

### Key Modeling Findings:
- **Regional Average Vulnerability:** The regional mean ADAVI score is **3.143 out of 5.000** (Moderate Vulnerability).
- **High to Extreme Vulnerability:** **2,501,292 acres (31.87% of total regional land)** fall into **High** (Class 4) or **Extreme** (Class 5) vulnerability tiers.
- **Kings County Hotspot:** Kings County exhibits the highest proportion of severe vulnerability in the region: **30.05% Extreme** and **37.17% High** (totaling **67.22%** in high/extreme stress), driven by intense summer thermal anomalies (>42°C) and deep Tulare Lake alluvial subsidence zones.
- **Priority Clusters:** **309 contiguous priority mitigation clusters** spanning **2.38 million acres** were delineated for focused groundwater recharge (MAR) and demand management.

---

## 2. Multi-Factor Scientific Framework & Weighting

The ADAVI model combines multi-spectral canopy moisture, thermal infrared land surface temperatures, vegetation vigor, alluvial plain micro-topography, and compaction modeling into a standardized multi-criteria weighted overlay framework:

$$\text{ADAVI} = \sum_{i=1}^{5} (w_i \times F_i)$$

$$\text{ADAVI} = 0.30 \times F_{\text{NDWI}} + 0.25 \times F_{\text{LST}} + 0.20 \times F_{\text{NDVI}} + 0.15 \times F_{\text{Slope}} + 0.10 \times F_{\text{Overdraft}}$$

| Factor | Parameter & Source | Weight ($w_i$) | Agronomic & Hydro-Geomorphic Rationale | 1–5 Reclassification Thresholds |
| :--- | :--- | :---: | :--- | :--- |
| **$F_1$** | **Crop Canopy Moisture (NDWI)**<br>*(Landsat 8/9 NIR/SWIR)* | **30%** | Measures water content in the crop canopy. Canopy desiccation directly signals irrigation shortfalls. | **1:** $>0.20$<br>**2:** $0.05 \text{ to } 0.20$<br>**3:** $-0.10 \text{ to } 0.05$<br>**4:** $-0.25 \text{ to } -0.10$<br>**5:** $\le -0.25$ |
| **$F_2$** | **Thermal Stress (LST °C)**<br>*(Landsat 8/9 TIRS)* | **25%** | Quantifies extreme evaporative demand and surface heating. Stressed, non-transpiring crops and fallowed soils reach peak thermal temperatures. | **1:** $\le 28^\circ\text{C}$<br>**2:** $28 \text{ to } 34^\circ\text{C}$<br>**3:** $34 \text{ to } 39^\circ\text{C}$<br>**4:** $39 \text{ to } 44^\circ\text{C}$<br>**5:** $>44^\circ\text{C}$ |
| **$F_3$** | **Crop Biomass / Vigor (NDVI)**<br>*(Landsat 8/9 Red/NIR)* | **20%** | Distinguishes productive irrigated orchards/row crops from fallowed, abandoned, or drought-stunted fields. | **1:** $>0.60$<br>**2:** $0.45 \text{ to } 0.60$<br>**3:** $0.30 \text{ to } 0.45$<br>**4:** $0.15 \text{ to } 0.30$<br>**5:** $\le 0.15$ |
| **$F_4$** | **Topographic Slope**<br>*(USGS 3DEP DEM)* | **15%** | Ultra-flat alluvial valley floors feature deep clay-silt layers (Corcoran Clay) prone to pore-collapse and compaction under groundwater pumping. | **1:** $>5.0^\circ$<br>**2:** $3.0 \text{ to } 5.0^\circ$<br>**3:** $1.8 \text{ to } 3.0^\circ$<br>**4:** $0.8 \text{ to } 1.8^\circ$<br>**5:** $\le 0.8^\circ$ |
| **$F_5$** | **Aquifer Overdraft Proxy**<br>*(Spatial Filtering & Geomorphology)* | **10%** | Delineates spatial corridors susceptible to inelastic aquifer compaction and continuous extraction cones of depression. | **1:** $<1.8$<br>**2:** $1.8 \text{ to } 2.6$<br>**3:** $2.6 \text{ to } 3.4$<br>**4:** $3.4 \text{ to } 4.2$<br>**5:** $\ge 4.2$ |

---

## 3. County-Level Vulnerability & Zonal Breakdown

The zonal statistics engine evaluated **Fresno**, **Tulare**, and **Kings** Counties independently:

| County | Total Area (Acres) | Mean ADAVI | Very Low (%) | Low (%) | Moderate (%) | High (%) | Extreme (%) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Kings County** | 885,173 | **3.828** | 0.00% | 6.11% | 26.66% | **37.17%** | **30.05%** |
| **Tulare County** | 3,090,288 | **3.032** | 18.82% | 26.81% | 30.73% | **15.59%** | **8.05%** |
| **Fresno County** | 3,798,971 | **3.083** | 23.49% | 16.98% | 28.74% | **19.14%** | **11.66%** |
| **Regional Total** | **7,848,462** | **3.143** | **18.79%** | **20.12%** | **29.22%** | **19.60%** | **12.27%** |

### County-Specific Observations:
1. **Kings County (Tulare Lake Bed Core):** Has by far the highest vulnerability concentration. With **595,082 acres (67.2%)** under High/Extreme vulnerability, the county experiences acute thermal extremes, almost no natural Sierra upland buffering within its borders, and high agricultural water demand.
2. **Fresno County (San Joaquin Basin):** Features **1,170,128 acres** in High/Extreme vulnerability, concentrated predominantly along the western agricultural belt (Westlands Water District / alluvial plains) where surface water allocations are often curtailed.
3. **Tulare County (Kaweah & Tule Basins):** Features **730,632 acres** in High/Extreme vulnerability, with severe stress pockets concentrated in the central-western valley floor and citrus/nut belts experiencing groundwater pumping stress.

---

## 4. Priority Groundwater Mitigation & Recharge Suitability

Using contiguous spatial clustering of Class 4 and 5 zones ($\ge 50$ hectares), **309 critical mitigation zones** were identified across **2,378,812 acres**:

```
[03_Drought_Aquifer_Model/Priority_Groundwater_Mitigation_Zones.geojson]
```

### Strategic Policy Recommendations:
1. **Managed Aquifer Recharge (MAR) & Flood-MAR Siting:** Channel excess winter atmospheric river flows onto high-permeability soils adjacent to identified high-vulnerability corridors to replenish depleted deep aquifers.
2. **Land Fallowing & Incentive Programs:** Target conservation easements and voluntary land repurposing programs (e.g., Multibenefit Land Repurposing Program) in Class 5 (Extreme Vulnerability) clusters to reduce pumping demand.
3. **Precision Crop Water Monitoring:** Deploy real-time satellite ET (e.g., OpenET) and thermal drone surveillance in High-Vulnerability zones to optimize variable-rate drip irrigation.
4. **SGMA Groundwater Allocation Caps:** Implement well-metering and tiered pumping fees focused on identified overdraft corridors to arrest inelastic ground subsidence.

---

## 5. Artifacts and Cartographic Deliverables

All datasets and cartographic products are cataloged in the project repository:
- **Processed Factor Rasters (`02_Processed_Factors/`):** Standardized 100m GeoTIFFs for Factors 1–5.
- **Model Output Rasters (`03_Drought_Aquifer_Model/`):** `ADAVI_Continuous_Index.tif` and `ADAVI_5Class_Vulnerability.tif`.
- **Zonal Statistics (`03_Drought_Aquifer_Model/`):** `County_Vulnerability_Statistics.csv` and `County_Vulnerability_Statistics.json`.
- **Mitigation Polygons (`03_Drought_Aquifer_Model/`):** `Priority_Groundwater_Mitigation_Zones.geojson`.
- **Publication Figures (`04_Final_Maps/`):**
  - `California_Central_Valley_Executive_4Panel.png` (300 DPI Executive Overview)
  - `California_Central_Valley_Cartographic_Map.png` (300 DPI Presentation Map with County Boundaries)
  - `Factor_Correlation_Matrix.png` (300 DPI Statistical Heatmap)
- **Interactive Notebook:** [`California_Central_Valley_Drought_Model.ipynb`](file:///c:/Users/USER/Documents/GIS/California_Central_Valley_Drought/California_Central_Valley_Drought_Model.ipynb).
