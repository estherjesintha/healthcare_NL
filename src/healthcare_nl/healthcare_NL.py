


# -------------------------
# 1. Imports
# -------------------------
import os
import requests
import pandas as pd
import geopandas as gpd
import json
import psycopg2
import osmnx as ox
from shapely.geometry import Point
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import numpy as np
from IPython.display import IFrame, display
import webbrowser
import folium
from pathlib import Path
import plotly.express as px
import getpass
import numpy as np
from shapely.geometry import box
# -------------------------
# 2. Database connection (interactive)
# -------------------------
# -------------------------
# Functions
# -------------------------

def download_wfs_geojson(url, layer_name, crs="EPSG:28992"):
    params = {
        "service": "WFS",
        "version": "1.1.0",
        "request": "GetFeature",
        "typeName": layer_name,
        "outputFormat": "application/json",
        "srsName": crs
    }
    r = requests.get(url, params=params)
    r.raise_for_status()
    return gpd.read_file(r.text)
def main():

        print("Enter PostgreSQL connection details:")

        db_host = input("Host (default: localhost): ").strip() or "localhost"
        db_port = input("Port (default: 5432): ").strip() or "5432"
        db_name = input("Database name: ").strip() or "postgres"
        db_user = input("Username: ").strip() or "postgres"
        db_password = getpass.getpass("Password: ") or "postgres"

        # Save into a config dictionary (optional but clean)
        db_config = {
            "host": db_host,
            "port": db_port,
            "database": db_name,
            "user": db_user,
            "password": db_password
        }

        # Try connecting
        try:
            conn = psycopg2.connect(**db_config)
            cur = conn.cursor()
            print("\n✅ Database connection successful!")

        except Exception as e:
            print("\n❌ Failed to connect to database:")
            print(e)
            raise

        # -------------------------
        # 3. Functions for WFS download
        # -------------------------
        def download_wfs_geojson(url, layer_name, crs="EPSG:28992"):
            """
            Fetches a WFS layer and returns it as a GeoDataFrame.
            """
            params = {
                "service": "WFS",
                "version": "1.1.0",
                "request": "GetFeature",
                "typeName": layer_name,
                "outputFormat": "application/json",
                "srsName": crs
            }
            r = requests.get(url, params=params)
            r.raise_for_status()
            return gpd.read_file(r.text)



        # -------------------------
        # 4. Download Administrative boundaries
        # -------------------------
        print("Downloading administrative boundaries (municipalities)...")
        muni = download_wfs_geojson(
            "https://service.pdok.nl/kadaster/bestuurlijkegebieden/wfs/v1_0?",
            "bg:Gemeentegebied"
        )
        muni.to_file("municipal_boundaries.geojson", driver="GeoJSON")
        
        # -------------------------
        # 5. Download Population data
        # -------------------------
        print("Downloading population data (CBS LAU)...")
        pop = download_wfs_geojson(
            "https://service.pdok.nl/cbs/pd/wfs/v1_0?",
            "pd:pd-nl-lau-2018"
        )
        pop.to_file("popdensitylau.geojson", driver="GeoJSON")

        # -------------------------
        # 4. Download Administrative boundaries (robust)
        # -------------------------
        print("Downloading administrative boundaries (municipalities)...")

        muni_path = "municipal_boundaries.geojson"

        if os.path.exists(muni_path):
            print("Loading cached municipal boundaries...")
            muni = gpd.read_file(muni_path)
        else:
            muni = download_wfs_geojson(
                "https://service.pdok.nl/kadaster/bestuurlijkegebieden/wfs/v1_0?",
                "bg:Gemeentegebied",
                crs="EPSG:4326"  # safer for WFS
            )
            if muni.empty:
                raise RuntimeError("Municipality WFS returned 0 features — service may be down or layer name changed.")
            muni = muni.to_crs(epsg=28992)
            muni.to_file(muni_path, driver="GeoJSON")

        print(f"Municipalities downloaded: {len(muni)}")
        print("CRS:", muni.crs)

        # -------------------------
        # 5. Download Population data (robust)
        # -------------------------
        print("\nDownloading population data (CBS LAU)...")

        pop_path = "popdensitylau.geojson"

        if os.path.exists(pop_path):
            print("Loading cached population data...")
            pop = gpd.read_file(pop_path)
        else:
            pop = download_wfs_geojson(
                "https://service.pdok.nl/cbs/pd/wfs/v1_0?",
                "pd:pd-nl-lau-2018",
                crs="EPSG:4326"
            )
            if pop.empty:
                raise RuntimeError("Population WFS returned 0 features — layer may be deprecated or temporarily unavailable.")
            pop = pop.to_crs(epsg=28992)
            pop.to_file(pop_path, driver="GeoJSON")

        print(f"Population polygons downloaded: {len(pop)}")
        print("CRS:", pop.crs)

        # -------------------------
        # Quick sanity checks
        # -------------------------
        print("\nSanity check:")
        print("Municipality columns:", list(muni.columns))
        print("Population columns:", list(pop.columns))

        # Check geometry validity
        print("Invalid municipality geometries:", (~muni.is_valid).sum())
        print("Invalid population geometries:", (~pop.is_valid).sum())

        # -------------------------
        # 6. Create tables in Postgres
        # -------------------------
        cur.execute("""
        CREATE TABLE IF NOT EXISTS public.municipalboundaries (
            id SERIAL PRIMARY KEY,
            identificatie text,
            naam text,
            code text,


            ligtinprovinciecode text,
            ligtinprovincienaam text,
            geom geometry(Geometry, 28992)
        )
        """)
        cur.execute("""
        CREATE TABLE IF NOT EXISTS public.popdensitylau (
            id SERIAL PRIMARY KEY,
            gml_id text,
            localid text,
            namespace text,
            versionid text,
            identifier text,
            identifierscheme text,
            language text,
            namestatus text,
            sourceofname text,
            pronunciation text,
            text text,
            script text,
            beginlifespanversion text,
            endlifespanversion text,
            timeposition text,
            referenceperiod_end text,
            leastdetailedscale integer,
            areavalue double precision,
            areavalue_uom text,
            landareavalue double precision,
            livableareavalue double precision,
            tesselation_localid text,
            tesselation_namespace text,
            pd_nl_lau_t_stat text,
            pd_nl_lau_t_measure text,
            pd_nl_lau_period_of_reference text,
            pd_nl_lau_t_obs_value integer,
            pd_nl_lau_d_stat text,
            pd_nl_lau_d_measure text,
            pd_nl_lau_d_obs_value double precision,
            pd_nl_lau_m_stat text,
            pd_nl_lau_m_measure text,
            pd_nl_lau_m_obs_value integer,
            pd_nl_lau_f_stat text,
            pd_nl_lau_f_measure text,
            pd_nl_lau_f_obs_value integer,
            geom geometry(Geometry, 28992)
        )
        """)
        conn.commit()

        

        # -------------------------
        # Ask user for province
        # -------------------------
        province_name = input("Enter province name (default: Zuid-Holland): ").strip()
        if province_name == "":
            province_name = "Zuid-Holland"
        print(f"Province selected: {province_name}")

        # -------------------------
        # Prepare province geometry
        # -------------------------
        province_munis = muni.to_crs(epsg=28992)
        province_munis = province_munis[
            province_munis["ligtInProvincieNaam"].str.lower() == province_name.lower()
        ].copy()

        if province_munis.empty:
            raise ValueError(f"No municipalities found for province '{province_name}'")

        # Dissolve municipalities into one polygon
        province_geom_m = province_munis.dissolve().geometry.iloc[0].buffer(0)

        # IMPORTANT: simplify in meters (reduces Overpass load)
        province_geom_m = province_geom_m.simplify(1000)  # 1 km tolerance

        # Convert to WGS84 for OSM query
        province_geom = gpd.GeoSeries([province_geom_m], crs="EPSG:28992").to_crs(epsg=4326).iloc[0]

        # -------------------------
        # Fetch healthcare facilities from OSM
        # -------------------------
        ox.settings.use_cache = True
        ox.settings.log_console = False

        health_tags = {"amenity": ["hospital", "clinic", "doctors"]}

        print("Querying OSM for healthcare facilities... (this may take ~30 seconds)")

        healthcare_gdf = ox.features_from_polygon(province_geom, health_tags)

        # Reproject to meters for analysis
        healthcare_gdf = healthcare_gdf.to_crs(epsg=28992)

        print(f"Number of healthcare facilities fetched in {province_name}: {len(healthcare_gdf)}")

        # Save output
        healthcare_gdf.to_file("healthcare_filtered.geojson", driver="GeoJSON")
        print("Saved to healthcare_filtered.geojson")

        # -------------------------
        # 7. Automatically derive extent and create fishnet (1 km)
        # -------------------------


        print("Creating 1 km fishnet over municipalities...")

        # Ensure CRS exists
        if muni.crs is None:
            raise ValueError("Municipal boundaries have no CRS defined.")

        # Reproject to metric CRS (RD New)
        muni_m = muni.to_crs(epsg=28992)

        # Get total bounds
        minx, miny, maxx, maxy = muni_m.total_bounds

        # Define grid resolution
        cell_size = 1000  # meters (1 km)

        # Generate grid cells
        xs = np.arange(minx, maxx, cell_size)
        ys = np.arange(miny, maxy, cell_size)

        grid_cells = []
        for x in xs:
            for y in ys:
                grid_cells.append(box(x, y, x + cell_size, y + cell_size))

        fishnet = gpd.GeoDataFrame(geometry=grid_cells, crs=muni_m.crs)

        # Clip grid to municipality boundaries
        fishnet = gpd.overlay(fishnet, muni_m, how="intersection")

        # Add unique ID
        fishnet["grid_id"] = range(len(fishnet))

        # Save fishnet
        fishnet.to_file("municipality_fishnet_1km.geojson", driver="GeoJSON")

        print(f"Fishnet created with {len(fishnet)} cells.")

        # ------------------------- #
        # Interactive map using Folium with title
        # ------------------------- #


        # Load files
        muni = gpd.read_file("municipal_boundaries.geojson")
        fishnet = gpd.read_file("municipality_fishnet_1km.geojson")

        # Convert to WGS84 for Folium
        muni_wgs = muni.to_crs(epsg=4326)
        fishnet_wgs = fishnet.to_crs(epsg=4326)

        # Compute map center
        center = muni_wgs.geometry.union_all().centroid
        center_latlon = [center.y, center.x]

        # Create map
        m = folium.Map(
            location=center_latlon,
            zoom_start=9,
            tiles="OpenStreetMap"
        )

        # Add municipality boundaries
        folium.GeoJson(
            muni_wgs,
            name="Municipal boundaries",
            style_function=lambda x: {
                "fillColor": "none",
                "color": "black",
                "weight": 2
            }
        ).add_to(m)

        # Add fishnet grid
        folium.GeoJson(
            fishnet_wgs,
            name="1 km Fishnet",
            style_function=lambda x: {
                "fillColor": "none",
                "color": "blue",
                "weight": 0.5
            }
        ).add_to(m)

        # Layer control
        folium.LayerControl().add_to(m)

        # ------------------------- #
        # Add title overlay
        # ------------------------- #
        title_html = """
        <div style="
        position: fixed;
        top: 10px;
        left: 50%;
        transform: translateX(-50%);
        z-index: 9999;
        background-color: white;
        padding: 10px 20px;
        border: 2px solid black;
        border-radius: 6px;
        font-size: 18px;
        font-weight: bold;
        box-shadow: 2px 2px 5px rgba(0,0,0,0.3);
        ">
        Municipality boundary with 1 km fishnet (EPSG:28992)
        </div>
        """

        m.get_root().html.add_child(folium.Element(title_html))

        # Save map
        map_file = Path("fishnet_map_with_title.html")
        m.save(map_file)

        print(f"Interactive map saved to: {map_file.resolve()}")
        map_file


        webbrowser.open(map_file.resolve().as_uri())

        # ------------------------- #
        #Compute centroids in EPSG:28992
        # ------------------------- #

        # Compute centroids
        fishnet["centroid"] = fishnet.geometry.centroid

        # Create centroid GeoDataFrame with only point geometry
        centroids = fishnet[["grid_id", "centroid"]].copy()
        centroids = centroids.set_geometry("centroid")
        centroids = centroids.rename_geometry("geometry")

        # Save centroids
        centroids.to_file("fishnet_centroids_1km_epsg28992.geojson", driver="GeoJSON")

        print("Centroids calculated and saved in EPSG:28992.")

        # ------------------------- #
        #1. Retrieve the number of healthcare facilities per province in the Netherlands.

        # ------------------------- #
        # Prepare province municipalities
        # ------------------------- #
        # province_name = "Zuid-Holland"

        province_munis = muni.to_crs(epsg=28992)
        province_munis = province_munis[
            province_munis["ligtInProvincieNaam"].str.lower() == province_name.lower()
        ].copy()

        if province_munis.empty:
            raise ValueError(f"No municipalities found for province '{province_name}'")

        # Dissolve province geometry
        province_geom_m = province_munis.dissolve().geometry.iloc[0].buffer(0)
        province_geom_m = province_geom_m.simplify(1000)

        province_geom = gpd.GeoSeries([province_geom_m], crs="EPSG:28992").to_crs(epsg=4326).iloc[0]

        # ------------------------- #
        # Fetch healthcare facilities
        # ------------------------- #
        ox.settings.use_cache = True
        ox.settings.log_console = False

        health_tags = {"amenity": ["hospital", "clinic", "doctors"]}

        print("Fetching healthcare facilities from OSM...")
        healthcare_zuid_holland = ox.features_from_polygon(province_geom, health_tags)
        healthcare_zuid_holland = healthcare_zuid_holland.to_crs(epsg=28992)

        print(f"Total healthcare facilities in {province_name}: {len(healthcare_zuid_holland)}")

        # ------------------------- #
        # Assign facilities to municipalities
        # ------------------------- #
        hc_muni = gpd.sjoin(
            healthcare_zuid_holland,
            province_munis[["naam", "geometry"]],
            how="left",
            predicate="within"
        )

        # ------------------------- #
        # Count facilities per municipality
        # ------------------------- #
        facility_count = (
            hc_muni.groupby("naam")
            .size()
            .reset_index(name="facility_count")
        )

        # ------------------------- #
        # Print top 5
        # ------------------------- #
        top5 = facility_count.sort_values("facility_count", ascending=False).head(5)
        print("\nQ 1. Retrieve the number of healthcare facilities per municipality in the province.")

        print("\nTop 5 municipalities by number of healthcare facilities:")
        for _, row in top5.iterrows():
            print(f"- {row['naam']}: {row['facility_count']} facilities")

        # ------------------------- #
        # Save outputs
        # ------------------------- #
        facility_count.to_csv("facility_count_by_municipality.csv", index=False)

        facility_count_gdf = province_munis.merge(facility_count, on="naam", how="left")
        facility_count_gdf["facility_count"] = facility_count_gdf["facility_count"].fillna(0)

        facility_count_gdf.to_file("facility_count_by_municipality.geojson", driver="GeoJSON")
        healthcare_zuid_holland.to_file("healthcare_zuid_holland.geojson", driver="GeoJSON")

        print("\nSaved files:")
        print("- healthcare_layer.geojson")
        print("- facility_count_by_municipality.csv")
        print("- facility_count_by_municipality.geojson")

        # ------------------------- #
        # 2. Calculate underserved area (>5 km to nearest healthcare) per municipality
        # -------------------------
        # Ensure CRS consistency
        # -------------------------
        fishnet = fishnet.to_crs(epsg=28992)
        healthcare_zuid_holland = healthcare_zuid_holland.to_crs(epsg=28992)
        province_munis = province_munis.to_crs(epsg=28992)

        # -------------------------
        # Build KD-tree from healthcare locations
        # -------------------------
        hc_points = np.array([
            (geom.x, geom.y)
            for geom in healthcare_zuid_holland.geometry
            if geom.geom_type == "Point"
        ])

        hc_tree = cKDTree(hc_points)

        # -------------------------
        # Compute distance from each fishnet centroid to nearest healthcare
        # -------------------------
        fishnet["centroid"] = fishnet.geometry.centroid

        centroid_coords = np.array([
            (p.x, p.y) for p in fishnet["centroid"]
        ])

        distances, _ = hc_tree.query(centroid_coords, k=1)
        fishnet["dist_to_hc"] = distances

        # -------------------------
        # Flag underserved cells (>5 km)
        # -------------------------
        fishnet["underserved"] = fishnet["dist_to_hc"] > 5000  # meters

        # -------------------------
        # Compute cell area
        # -------------------------
        fishnet["cell_area"] = fishnet.geometry.area

        # -------------------------
        # Spatial join fishnet → municipality
        # -------------------------
        fishnet_muni = fishnet.sjoin(
            province_munis[["naam", "geometry"]],
            how="left",
            predicate="within"
        )

        # -------------------------
        # Aggregate underserved area per municipality
        # -------------------------
        underserved_area = (
            fishnet_muni[fishnet_muni["underserved"]]
            .groupby("naam_right")["cell_area"]
            .sum()
            .reset_index()
            .rename(columns={
                "naam_right": "naam",
                "cell_area": "underserved_area_m2"
            })
        )

        # -------------------------
        # Merge back and compute km²
        # -------------------------
        result = province_munis.merge(underserved_area, on="naam", how="left")
        result["underserved_area_m2"] = result["underserved_area_m2"].fillna(0)
        result["underserved_area_km2"] = result["underserved_area_m2"] / 1e6

        # -------------------------
        # Output (Top 5 by underserved area)
        # -------------------------
        print(
            "\nQ 2. Top 5 municipalities by area beyond the recommended 5 km walking distance to healthcare:"
        )

        top5_underserved = (
            result[["naam", "underserved_area_km2"]]
            .sort_values("underserved_area_km2", ascending=False)
            .head(5)
        )

        print(top5_underserved.to_string(index=False))

        # Save full layer
        result.to_file("zuid_holland_underserved_area.geojson", driver="GeoJSON")
        print("\nSaved to zuid_holland_underserved_area.geojson")



        # -------------------------
        # Safety checks
        # -------------------------
        if "underserved_area_km2" not in result.columns:
            raise ValueError("Column 'underserved_area_km2' not found in result.")

        # Replace NaNs just for visualization
        result_plot = result.copy()
        result_plot["underserved_area_km2"] = result_plot["underserved_area_km2"].fillna(0)

        # -------------------------
        # Convert to WGS84
        # -------------------------
        result_wgs = result_plot.to_crs(epsg=4326)
        hc_wgs = healthcare_zuid_holland.to_crs(epsg=4326)

        # -------------------------
        # Center map
        # -------------------------
        center = result_wgs.geometry.union_all().centroid
        center_latlon = [center.y, center.x]

        # -------------------------
        # Create map
        # -------------------------
        m = folium.Map(location=center_latlon, zoom_start=9, tiles="OpenStreetMap")

        # -------------------------
        # Color scale function
        # -------------------------
        def color_scale(val):
            if val == 0:
                return "#f0f0f0"  # no underserved
            elif val < 2:
                return "#fee8c8"
            elif val < 5:
                return "#fdbb84"
            else:
                return "#e34a33"

        # -------------------------
        # Add municipalities
        # -------------------------
        folium.GeoJson(
            result_wgs,
            name="Underserved Area",
            style_function=lambda feat: {
                "fillColor": color_scale(feat["properties"].get("underserved_area_km2", 0)),
                "color": "black",
                "weight": 0.8,
                "fillOpacity": 0.7,
            },
            tooltip=folium.GeoJsonTooltip(
                fields=["naam", "underserved_area_km2"],
                aliases=["Municipality:", "Underserved area (km²):"],
                localize=True,
                sticky=True
            )
        ).add_to(m)

        # -------------------------
        # Add healthcare points
        # -------------------------
        for _, row in hc_wgs.iterrows():
            geom = row.geometry
            if geom is None:
                continue

            if geom.geom_type == "Point":
                latlon = [geom.y, geom.x]
            else:
                c = geom.centroid
                latlon = [c.y, c.x]

            folium.CircleMarker(
                location=latlon,
                radius=3,
                color="blue",
                fill=True,
                fill_opacity=0.7,
            ).add_to(m)

        # -------------------------
        # Add title
        # -------------------------
        title_html = """
        <h3 align="center" style="font-size:20px">
        <b>Underserved Healthcare Areas </b>
        </h3>
        """
        m.get_root().html.add_child(folium.Element(title_html))

        # -------------------------
        # Add legend
        # -------------------------
        legend_html = """
        <div style="
        position: fixed;
        bottom: 40px;
        left: 40px;
        width: 260px;
        background-color: white;
        border:2px solid grey;
        z-index:9999;
        font-size:14px;
        padding:10px;">
        <b>Legend</b><br>
        <span style="background-color:#f0f0f0;width:12px;height:12px;display:inline-block;"></span> 0 km² underserved<br>
        <span style="background-color:#fee8c8;width:12px;height:12px;display:inline-block;"></span> &lt; 2 km² underserved<br>
        <span style="background-color:#fdbb84;width:12px;height:12px;display:inline-block;"></span> 2–5 km² underserved<br>
        <span style="background-color:#e34a33;width:12px;height:12px;display:inline-block;"></span> &gt; 5 km² underserved<br>
        <span style="color:blue;">●</span> Healthcare facility
        </div>
        """
        m.get_root().html.add_child(folium.Element(legend_html))

        # -------------------------
        # Save & open
        # -------------------------
        map_file = Path("underserved_healthcare_map.html")
        m.save(map_file)

        print(f"Map saved to: {map_file.resolve()}")
        webbrowser.open(map_file.resolve().as_uri())

        # ------------------------- #
        # 3. Calculate underserved ratio per municipality
        # ------------------------- #

        # Compute total municipality area
        result["total_area_m2"] = result.geometry.area

        # Compute underserved ratio = underserved area / total area
        result["underserved_ratio"] = result["underserved_area_m2"] / result["total_area_m2"]

        # Clean invalid values
        result["underserved_ratio"] = result["underserved_ratio"].replace([np.inf, -np.inf], np.nan).fillna(0)

        # ------------------------- #
        # Sort and print Top 5
        # ------------------------- #
        print("\nQ 3. Calculate the ratio of underserved area to total municipality area.")
        top5_ratio = result.sort_values("underserved_ratio", ascending=False).head(5)

        print("\nTop 5 municipalities by underserved area ratio:")
        for _, row in top5_ratio.iterrows():
            print(f"- {row['naam']}: {row['underserved_ratio']:.3f}")

        # -------------------------
        #4. Report the municipalities with the lowest and highest underserved area.
        # Municipality with lowest underserved area
        # -------------------------
        lowest = result.loc[result["underserved_area_m2"].idxmin()]

        # -------------------------
        # Municipality with highest underserved area
        # -------------------------
        print("\nQ 4. Report the municipalities with the lowest and highest underserved area.")
        highest = result.loc[result["underserved_area_m2"].idxmax()]

        print("\nMunicipality with LOWEST underserved area:")
        print(f" - {lowest['naam']}: {lowest['underserved_area_m2'] / 1e6:.2f} km²")

        print("\nMunicipality with HIGHEST underserved area:")
        print(f" - {highest['naam']}: {highest['underserved_area_m2'] / 1e6:.2f} km²")

        # ------------------------- #
        # 5. Identify municipalities with high population density but low healthcare facility density
        # (robust, safe to re-run)
        # ------------------------- #
        pop["PD_NL_LAU_T_OBS_VALUE"] = pd.to_numeric(pop["PD_NL_LAU_T_OBS_VALUE"], errors="coerce")
        pop["landAreaValue"] = pd.to_numeric(pop["landAreaValue"], errors="coerce")

        print("Starting robust high-risk municipality analysis...")

        # -------------------------
        # Safety checks
        # -------------------------
        required = ["province_munis", "healthcare_zuid_holland", "pop"]
        for r in required:
            if r not in locals():
                raise NameError(f"Missing required object: {r}")

        province_munis = province_munis.to_crs(epsg=28992)
        healthcare_zuid_holland = healthcare_zuid_holland.to_crs(epsg=28992)
        pop = pop.to_crs(epsg=28992)

        # -------------------------
        # Auto-detect population & area columns
        # -------------------------
        pop_col = next((c for c in pop.columns if "obs" in c.lower()), None)
        area_col = next((c for c in pop.columns if "area" in c.lower()), None)

        if pop_col is None or area_col is None:
            raise ValueError("Could not detect population or area column automatically.")

        print(f"Using population column: {pop_col}")
        print(f"Using area column: {area_col}")

        # -------------------------
        # Clean population table
        # -------------------------
        pop[pop_col] = pd.to_numeric(pop[pop_col], errors="coerce")
        pop[area_col] = pd.to_numeric(pop[area_col], errors="coerce")

        pop = pop.dropna(subset=[pop_col, area_col])
        pop = pop[pop[area_col] > 0]

        if pop.empty:
            raise ValueError("Population table empty after cleaning.")

        # -------------------------
        # Spatial join population → municipalities
        # -------------------------
        pop_muni = gpd.sjoin(
            pop,
            province_munis[["naam", "geometry"]],
            how="inner",
            predicate="intersects"
        )

        print(f"Population polygons assigned to municipalities: {len(pop_muni)}")

        # -------------------------
        # Compute population density per municipality
        # -------------------------
        pop_density = (
            pop_muni.groupby("naam")[[pop_col, area_col]]
            .sum()
            .reset_index()
        )

        pop_density["pop_density"] = pop_density[pop_col] / pop_density[area_col]

        # -------------------------
        # Build result base table
        # -------------------------
        result = province_munis.copy()

        # Merge pop density
        result = result.merge(pop_density[["naam", "pop_density"]], on="naam", how="left")
        result["pop_density"] = result["pop_density"].replace([np.inf, -np.inf], np.nan)

        # -------------------------
        # Healthcare density
        # -------------------------
        hc_muni = gpd.sjoin(
            healthcare_zuid_holland,
            province_munis[["naam", "geometry"]],
            how="inner",
            predicate="within"
        )

        hc_count = hc_muni.groupby("naam").size().reset_index(name="facility_count")

        result = result.merge(hc_count, on="naam", how="left")
        result["facility_count"] = result["facility_count"].fillna(0)

        result["area_km2"] = result.geometry.area / 1e6
        result["facility_density"] = result["facility_count"] / result["area_km2"]
        result["facility_density"] = result["facility_density"].replace([np.inf, -np.inf], np.nan)

        # -------------------------
        # Thresholds
        # -------------------------
        if result["pop_density"].dropna().empty:
            raise ValueError("No valid population density values found.")

        pop_thresh = result["pop_density"].dropna().quantile(0.75)
        facility_thresh = result["facility_density"].dropna().quantile(0.25)

        # -------------------------
        # High-risk classification
        # -------------------------
        result["high_risk"] = (
            (result["pop_density"] >= pop_thresh) &
            (result["facility_density"] <= facility_thresh)
        )

        high_risk = result[result["high_risk"]]

        # -------------------------
        # Output
        # -------------------------
        print("\nQ 5. Identify municipalities healthcare facility imbalance, i.e., high population density but low healthcare facility density.")
        print(f"\nPopulation density threshold (75th percentile): {pop_thresh:.4f}")
        print(f"Facility density threshold (25th percentile): {facility_thresh:.4f}")

        print("\nHigh-risk municipalities (high demand, low supply):")
        if high_risk.empty:
            print("None identified.")
        else:
            for _, r in high_risk.iterrows():
                print(f"- {r['naam']}: PopDensity={r['pop_density']:.4f}, FacilityDensity={r['facility_density']:.4f}")

        print("\nAnalysis complete.")

        # -------------------------
        # map visualization
        # -------------------------
        # -------------------------
        # Safety checks
        # -------------------------
        required = ["result", "healthcare_zuid_holland"]
        for name in required:
            if name not in locals():
                raise NameError(f"Required object '{name}' not found. Run previous steps first.")

        # -------------------------
        # Convert to WGS84 for folium
        # -------------------------
        result_wgs = result.to_crs(epsg=4326)
        hc_wgs = healthcare_zuid_holland.to_crs(epsg=4326)

        # -------------------------
        # Get map center
        # -------------------------
        center = result_wgs.geometry.union_all().centroid
        center_latlon = [center.y, center.x]

        # -------------------------
        # Create base map
        # -------------------------
        m = folium.Map(
            location=center_latlon,
            zoom_start=9,
            tiles="OpenStreetMap"
        )

        m = folium.Map(
            location=center_latlon,
            zoom_start=9,
            tiles="OpenStreetMap"
        )

        # Add title
        title_html = """
        <h3 align="center" style="
            font-size:22px;
            font-weight:bold;
            margin-top:10px;
        ">
        High-Risk Healthcare Access Map
        </h3>
        """
        m.get_root().html.add_child(folium.Element(title_html))


        # -------------------------
        # Add municipalities
        # -------------------------
        for _, row in result_wgs.iterrows():
            color = "red" if row.get("high_risk", False) else "#cccccc"

            folium.GeoJson(
                row.geometry,
                style_function=lambda x, col=color: {
                    "fillColor": col,
                    "color": "black",
                    "weight": 1,
                    "fillOpacity": 0.6 if col == "red" else 0.3
                },
                tooltip=f"{row['naam']}<br>Pop density: {row.get('pop_density', 'NA'):.4f}<br>Facility density: {row.get('facility_density', 'NA'):.4f}"
            ).add_to(m)

        # -------------------------
        # Add healthcare facilities safely
        # -------------------------
        for _, row in hc_wgs.iterrows():
            geom = row.geometry

            # Convert polygon → centroid if needed
            if geom.geom_type == "Point":
                latlon = [geom.y, geom.x]
            else:
                c = geom.centroid
                latlon = [c.y, c.x]

            folium.CircleMarker(
                location=latlon,
                radius=3,
                color="blue",
                fill=True,
                fill_opacity=0.7,
                popup=row.get("name", "Healthcare facility")
            ).add_to(m)

        # -------------------------
        # Add legend (HTML)
        # -------------------------
        legend_html = """
        <div style="
        position: fixed;
        bottom: 50px;
        left: 50px;
        width: 220px;
        height: 120px;
        background-color: white;
        border:2px solid grey;
        z-index:9999;
        font-size:14px;
        padding:10px;">
        <b>Legend</b><br>
        <span style="background-color:red;width:12px;height:12px;display:inline-block;"></span> High-risk municipality<br>
        <span style="background-color:#cccccc;width:12px;height:12px;display:inline-block;"></span> Other municipality<br>
        <span style="color:blue;">●</span> Healthcare facility
        </div>
        """
        m.get_root().html.add_child(folium.Element(legend_html))

        # -------------------------
        # Save and display
        # -------------------------
        map_file = Path("high_risk_healthcare_map.html")
        m.save(map_file)

        print(f"Map saved to: {map_file.resolve()}")
        map_file



        # Save map
        map_file = Path("high_risk_healthcare_map.html")
        m.save(map_file)

        # # Display inline inside Jupyter
        # display(IFrame(src=str(map_file), width="100%", height="600"))

        # Open in default browser (Windows)
        webbrowser.open(map_file.resolve().as_uri())

        # -------------------------
        # Scatter plot visualization (force open in browser)
        # -------------------------


        # Build dataframe for plotting
        plot_df = result[["pop_density", "facility_density", "high_risk", "naam"]].copy()
        plot_df["Risk"] = plot_df["high_risk"].map({True: "High risk", False: "Normal"})

        # Create interactive scatter plot
        fig = px.scatter(
            plot_df,
            x="pop_density",
            y="facility_density",
            color="Risk",
            hover_name="naam",
            title="Demand vs Supply by Municipality",
            labels={
                "pop_density": "Population Density",
                "facility_density": "Healthcare Facility Density"
            }
        )

        # Add threshold lines
        fig.add_vline(x=pop_thresh, line_dash="dash", line_color="red", annotation_text="High demand threshold")
        fig.add_hline(y=facility_thresh, line_dash="dash", line_color="blue", annotation_text="Low supply threshold")

        # Save as HTML
        html_file = Path("demand_vs_supply_scatter.html")
        fig.write_html(html_file)

        # Open in browser
        webbrowser.open(html_file.resolve().as_uri())



if __name__ == "__main__":
    main()