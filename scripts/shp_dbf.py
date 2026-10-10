"""Minimal pure-Python .shp / .dbf readers (Polygon shape type only).

No third-party geospatial deps available in this sandbox (pypi/apt egress is
blocked by org policy), so this hand-rolls just enough of the ESRI Shapefile
and dBASE formats to get county polygons + attributes out of the Census
cartographic boundary file.
"""
import struct


def read_dbf(path):
    """Return (field_names, list_of_row_dicts) for a .dbf file."""
    with open(path, "rb") as f:
        data = f.read()
    n_records = struct.unpack("<I", data[4:8])[0]
    header_size = struct.unpack("<H", data[8:10])[0]
    record_size = struct.unpack("<H", data[10:12])[0]

    fields = []
    offset = 32
    while data[offset] != 0x0D:
        raw = data[offset:offset + 32]
        name = raw[0:11].split(b"\x00")[0].decode("ascii")
        length = raw[16]
        fields.append((name, length))
        offset += 32

    rows = []
    rec_start = header_size
    for i in range(n_records):
        rec = data[rec_start:rec_start + record_size]
        rec_start += record_size
        if not rec or rec[0:1] == b"*":
            continue
        pos = 1
        row = {}
        for name, length in fields:
            raw = rec[pos:pos + length]
            row[name] = raw.decode("ascii", errors="replace").strip()
            pos += length
        rows.append(row)
    return [f[0] for f in fields], rows


def read_shp_polygons(path):
    """Yield one list-of-rings per feature, in file order. Each ring is a
    list of (x, y) tuples. Rings belonging to the same feature (multipart
    polygons, e.g. islands or donut holes) are all in the same list."""
    with open(path, "rb") as f:
        data = f.read()

    shapetype = struct.unpack("<i", data[32:36])[0]
    assert shapetype == 5, f"expected Polygon (5), got {shapetype}"

    pos = 100
    n = len(data)
    features = []
    while pos < n:
        rec_num = struct.unpack(">i", data[pos:pos + 4])[0]
        content_len_words = struct.unpack(">i", data[pos + 4:pos + 8])[0]
        content_len_bytes = content_len_words * 2
        content = data[pos + 8:pos + 8 + content_len_bytes]
        pos += 8 + content_len_bytes

        rec_shapetype = struct.unpack("<i", content[0:4])[0]
        if rec_shapetype == 0:
            features.append([])
            continue
        num_parts = struct.unpack("<i", content[36:40])[0]
        num_points = struct.unpack("<i", content[40:44])[0]
        parts_start = 44
        parts = struct.unpack(f"<{num_parts}i", content[parts_start:parts_start + 4 * num_parts])
        points_start = parts_start + 4 * num_parts
        pts = struct.unpack(f"<{2 * num_points}d", content[points_start:points_start + 16 * num_points])
        xy = list(zip(pts[0::2], pts[1::2]))

        rings = []
        for i, start_idx in enumerate(parts):
            end_idx = parts[i + 1] if i + 1 < len(parts) else num_points
            rings.append(xy[start_idx:end_idx])
        features.append(rings)
    return features


def point_in_rings(px, py, rings):
    """Even-odd rule across all rings of a (possibly multipart, possibly
    donut-holed) polygon feature. Standard ray-casting per ring, XORed."""
    inside = False
    for ring in rings:
        n = len(ring)
        if n < 3:
            continue
        j = n - 1
        for i in range(n):
            xi, yi = ring[i]
            xj, yj = ring[j]
            if ((yi > py) != (yj > py)) and (
                px < (xj - xi) * (py - yi) / (yj - yi + 1e-300) + xi
            ):
                inside = not inside
            j = i
    return inside


def ring_bbox(rings):
    xs = [p[0] for ring in rings for p in ring]
    ys = [p[1] for ring in rings for p in ring]
    if not xs:
        return None
    return (min(xs), min(ys), max(xs), max(ys))
