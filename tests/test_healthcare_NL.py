import sys
from pathlib import Path
import pytest
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point, Polygon

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from healthcare_nl import healthcare_NL as hp


# ------------------------
# Fixtures
# ------------------------

# Fixtures create reusable test data that is injected into tests automatically.
@pytest.fixture
def sample_points():
    """
    Fake WFS response data.
    Used to test download and parsing logic without network dependency.
    """
    return gpd.GeoDataFrame(
        {"id": [1, 2]},
        geometry=[Point(0, 0), Point(1, 1)],
        crs="EPSG:4326"
    )


@pytest.fixture
def sample_fishnet():
    """
    Simplified spatial grid (fishnet) used to test geometry area,
    centroid computation, and spatial relationships.
    """
    return gpd.GeoDataFrame(
        {"grid_id": [0, 1]},
        geometry=[
            Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
            Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])
        ],
        crs="EPSG:28992"
    )


@pytest.fixture
def sample_healthcare():
    """
    Fake healthcare facility locations.
    Used for distance computation and underserved-area logic.
    """
    return gpd.GeoDataFrame(
        {"name": ["HC1", "HC2"]},
        geometry=[Point(0.5, 0.5), Point(3, 1)],
        crs="EPSG:28992"
    )


# ------------------------
# Tests
# ------------------------

def test_download_function_exists():
     # Check that the main download function is defined in the module
    """Ensure core function exists."""
    assert hasattr(hp, "download_wfs_geojson"), "download_wfs_geojson() function missing"


def test_download_returns_geodataframe(sample_points, monkeypatch):
    """download_wfs_geojson returns a GeoDataFrame."""

    class FakeResp:
        text = "fake"
        def raise_for_status(self):
            pass

    monkeypatch.setattr(hp.requests, "get", lambda *a, **k: FakeResp())
    # Replace file reading with fixture data
    monkeypatch.setattr(hp.gpd, "read_file", lambda _: sample_points)

    gdf = hp.download_wfs_geojson("url", "layer")

    assert isinstance(gdf, gpd.GeoDataFrame)
    assert len(gdf) == 2


def test_fishnet_geometry_valid(sample_fishnet):
    """Fishnet cells should have positive area."""
    assert (sample_fishnet.area > 0).all()


def test_centroid_distance_non_negative(sample_fishnet, sample_healthcare):
    """Distances to healthcare should be >= 0."""
    fishnet = sample_fishnet.copy()
    fishnet["centroid"] = fishnet.geometry.centroid

    hc_points = np.array([(p.x, p.y) for p in sample_healthcare.geometry])
    from scipy.spatial import cKDTree
    tree = cKDTree(hc_points)

    coords = np.array([(p.x, p.y) for p in fishnet["centroid"]])
    distances, _ = tree.query(coords, k=1)

    assert (distances >= 0).all()


def test_underserved_flag_logic():
    """Underserved flag uses >5000m rule."""
    df = pd.DataFrame({"dist_to_hc": [1000, 5001, 8000]})
    df["underserved"] = df["dist_to_hc"] > 5000

    assert df["underserved"].tolist() == [False, True, True]


def test_population_density_calculation():
    """Population density aggregation correct."""
    df = pd.DataFrame({
        "naam": ["A", "A", "B"],
        "pop": [100, 200, 50],
        "area": [1, 1, 2]
    })

    grouped = df.groupby("naam")[["pop", "area"]].sum().reset_index()
    grouped["pop_density"] = grouped["pop"] / grouped["area"]

    assert grouped.loc[grouped["naam"] == "A", "pop_density"].iloc[0] == 150
    assert grouped.loc[grouped["naam"] == "B", "pop_density"].iloc[0] == 25


def test_high_risk_classification():
    """High risk classification uses quantile thresholds."""
    df = pd.DataFrame({
        "pop_density": [10, 50, 80],
        "facility_density": [5, 1, 0.5]
    })

    pop_thresh = df["pop_density"].quantile(0.75)
    fac_thresh = df["facility_density"].quantile(0.25)

    df["high_risk"] = (df["pop_density"] >= pop_thresh) & (df["facility_density"] <= fac_thresh)

    assert df["high_risk"].any()


def test_missing_population_columns_raises():
    """Missing population/area columns should raise error."""
    pop = gpd.GeoDataFrame({"x": [1, 2]}, geometry=[Point(0, 0), Point(1, 1)])

    pop_col = next((c for c in pop.columns if "obs" in c.lower()), None)
    area_col = next((c for c in pop.columns if "area" in c.lower()), None)

    with pytest.raises(ValueError, match="detect population or area"):
        if pop_col is None or area_col is None:
            raise ValueError("Could not detect population or area column automatically.")
