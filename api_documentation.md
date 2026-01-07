# API Documentation

This file documents the APIs and web services used in the project. It provides endpoints, layers, and purpose for each dataset.

---

## 1. Municipal Boundaries (Administrative Areas)

- **API Type:** RESTful WFS (Web Feature Service)
- **Endpoint:** `https://service.pdok.nl/kadaster/bestuurlijkegebieden/wfs/v1_0?`
- **Layer Name:** `bg:Gemeentegebied`
- **Description:** Provides boundaries of all municipalities in the Netherlands. Includes attributes such as:
  - `naam`: Municipality name
  - `ligtInProvincieNaam`: Province the municipality belongs to
  - Geometry in polygon format
- **Purpose in Project:**
  - Define municipal boundaries for analysis and map visualizations
  - Used in creating fishnet grids for distance calculations to healthcare facilities
- **CRS:** Retrieved in EPSG:4326, converted to EPSG:28992 (Amersfoort / RD New)

---

## 2. Population per Municipality

- **API Type:** RESTful WFS (Web Feature Service)
- **Endpoint:** `https://service.pdok.nl/cbs/pd/wfs/v1_0?`
- **Layer Name:** `pd:pd-nl-lau-2018`
- **Description:** Provides population and area data for municipalities. Key attributes:
  - Population counts (`PD_NL_LAU_T_OBS_VALUE`)
  - Land area (`landAreaValue`)
  - Geometry in polygons corresponding to municipal areas
- **Purpose in Project:**
  - Calculate population density per municipality
  - Identify high-risk municipalities with high population but low healthcare facility coverage
- **CRS:** Retrieved in EPSG:4326, converted to EPSG:28992

---

## 3. OpenStreetMap Healthcare Facilities

- **API Type:** Overpass / OSMnx Python library
- **Data Tags:** `amenity=hospital`, `amenity=clinic`, `amenity=doctors`
- **Description:** Provides point locations of healthcare facilities across selected provinces.
- **Purpose in Project:**
  - Calculate facility density per municipality
  - Compute distances to nearest healthcare for underserved area analysis
- **CRS:** Retrieved in WGS84 (EPSG:4326), converted to EPSG:28992 for analysis

---

## Notes

1. All spatial data is rendered in EPSG:28992 for consistent distance calculations and mapping.
2. Fishnet grids (1 km x 1 km) are generated for accurate distance-based underserved area estimation.
3. The combination of population, administrative boundaries, and healthcare facility data enables:
   - Identifying underserved municipalities
   - Visualizing high-risk areas on interactive maps
   - Calculating demand vs supply ratios

