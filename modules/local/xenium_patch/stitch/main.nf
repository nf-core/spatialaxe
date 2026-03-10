/*
 * XENIUM_PATCH_STITCH: Stitch per-patch segmentation results into unified output.
 *
 * Uses sopa's solve_conflicts() to resolve overlapping cells at patch boundaries.
 *
 * Input:
 *   - meta: Sample metadata map
 *   - patches: Directory containing patch subdirectories and patch_grid.json
 *
 * Output:
 *   - xr_polygons_transcript: Stitched cell polygons and transcript metadata
 *   - versions: Software versions
 */
process XENIUM_PATCH_STITCH {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9c8f3a2de4e2aa94500011f7d7d09276e9b6f2d79ee8737c9098fe22d4649bc/data' :
        'community.wave.seqera.io/library/sopa_procps-ng_pyarrow:c9ce8cd2ede79d72' }"

    input:
    tuple val(meta), path(patches)

    output:
    tuple val(meta),
        path("output/xr-cell-polygons.geojson"),
        path("output/xr-transcript-metadata.csv")  , emit: xr_polygons_transcript
    path("versions.yml")                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    stitch_transcripts.py \\
        --patches ${patches} \\
        --output output \\
        --min-transcripts-per-cell ${params.baysor_tiling_min_transcripts_per_cell}

    # Post-process: ensure all GeoJSON geometries are Polygon.
    # make_valid() and solve_conflicts() can produce MultiPolygon,
    # MultiLineString, or GeometryCollection — XeniumRanger rejects these.
    # Dropped cells must also be removed from the transcript CSV.
    python3 -c "
import csv, json, shapely
from shapely.geometry import mapping, shape

geojson_path = 'output/xr-cell-polygons.geojson'
csv_path = 'output/xr-transcript-metadata.csv'

with open(geojson_path) as f:
    data = json.load(f)

clean = []
dropped_cells = set()
for feat in data['features']:
    geom = shape(feat['geometry'])
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    poly = None
    if geom.geom_type == 'Polygon':
        poly = geom
    elif geom.geom_type == 'MultiPolygon':
        poly = max(geom.geoms, key=lambda g: g.area)
    elif geom.geom_type == 'GeometryCollection':
        polys = [g for g in geom.geoms if g.geom_type == 'Polygon']
        if polys:
            poly = max(polys, key=lambda g: g.area)
    if poly is not None and not poly.is_empty:
        feat['geometry'] = mapping(poly)
        clean.append(feat)
    else:
        cell_id = feat.get('id') or feat.get('properties', {}).get('cell_id', '')
        dropped_cells.add(str(cell_id))

print(f'GeoJSON: {len(clean)} kept, {len(dropped_cells)} dropped: {dropped_cells}')
data['features'] = clean
with open(geojson_path, 'w') as f:
    json.dump(data, f)

if dropped_cells:
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    reassigned = 0
    for row in rows:
        if row['cell'] in dropped_cells:
            row['cell'] = ''
            row['is_noise'] = '1'
            reassigned += 1
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f'CSV: {reassigned} transcripts reassigned to UNASSIGNED')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        sopa: \$(python3 -c "import sopa; print(sopa.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p output
    echo '{"type":"FeatureCollection","features":[]}' > output/xr-cell-polygons.geojson
    echo 'transcript_id,x,y,z,gene,cell,is_noise' > output/xr-transcript-metadata.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.12.0"
        sopa: "0.1.0"
    END_VERSIONS
    """
}
