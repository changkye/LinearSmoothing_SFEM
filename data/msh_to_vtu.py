#!/usr/bin/env python3

from __future__ import annotations

import argparse
import pathlib
from typing import Dict, List, Tuple


def gmsh_nodes_per_element(element_type: int) -> int:
    mapping = {
        1: 2,   # line
        2: 3,   # triangle
        8: 3,   # second-order line
        9: 6,   # second-order triangle
        15: 1,  # point
    }
    if element_type not in mapping:
        raise ValueError(f"Unsupported Gmsh element type: {element_type}")
    return mapping[element_type]


def parse_msh_ascii_v4(file_path: pathlib.Path) -> Tuple[List[Tuple[float, float, float]], List[List[int]]]:
    tokens = file_path.read_text(encoding="utf-8").split()
    idx = 0

    def require(expected: str) -> None:
        nonlocal idx
        if idx >= len(tokens) or tokens[idx] != expected:
            raise ValueError(f"Malformed msh file {file_path}: expected {expected}")
        idx += 1

    require("$MeshFormat")
    version = float(tokens[idx])
    file_type = int(tokens[idx + 1])
    idx += 3
    if version < 4.0:
        raise ValueError(f"Unsupported msh version in {file_path}: {version}")
    if file_type != 0:
        raise ValueError(f"Only ASCII msh is supported: {file_path}")
    require("$EndMeshFormat")

    node_tags_to_index: Dict[int, int] = {}
    nodes: List[Tuple[float, float, float]] = []
    elements: List[List[int]] = []

    while idx < len(tokens):
        section = tokens[idx]
        idx += 1

        if section == "$Nodes":
            num_entity_blocks = int(tokens[idx])
            num_nodes = int(tokens[idx + 1])
            idx += 4

            for _ in range(num_entity_blocks):
                entity_dim = int(tokens[idx])
                _entity_tag = int(tokens[idx + 1])
                parametric = int(tokens[idx + 2])
                num_nodes_in_block = int(tokens[idx + 3])
                idx += 4

                block_tags = [int(tokens[idx + i]) for i in range(num_nodes_in_block)]
                idx += num_nodes_in_block

                for tag in block_tags:
                    x = float(tokens[idx])
                    y = float(tokens[idx + 1])
                    z = float(tokens[idx + 2])
                    idx += 3
                    for _ in range(parametric * entity_dim):
                        idx += 1
                    node_tags_to_index[tag] = len(nodes)
                    nodes.append((x, y, z))

            require("$EndNodes")
            if len(nodes) != num_nodes:
                raise ValueError(f"Node count mismatch while reading {file_path}")
            continue

        if section == "$Elements":
            num_entity_blocks = int(tokens[idx])
            _num_elements = int(tokens[idx + 1])
            idx += 4

            for _ in range(num_entity_blocks):
                _entity_dim = int(tokens[idx])
                _entity_tag = int(tokens[idx + 1])
                element_type = int(tokens[idx + 2])
                num_elements_in_block = int(tokens[idx + 3])
                idx += 4

                nodes_per_element = gmsh_nodes_per_element(element_type)
                for _ in range(num_elements_in_block):
                    _element_tag = int(tokens[idx])
                    idx += 1
                    raw_nodes = [int(tokens[idx + i]) for i in range(nodes_per_element)]
                    idx += nodes_per_element
                    if element_type == 9:
                        elements.append([node_tags_to_index[tag] for tag in raw_nodes])

            require("$EndElements")
            continue

        if section.startswith("$"):
            end_section = "$End" + section[1:]
            while idx < len(tokens) and tokens[idx] != end_section:
                idx += 1
            require(end_section)
            continue

        raise ValueError(f"Unexpected token while parsing {file_path}: {section}")

    if not nodes:
        raise ValueError(f"No nodes found in {file_path}")
    if not elements:
        raise ValueError(f"No T6 triangle elements found in {file_path}")
    return nodes, elements


def write_vtu(file_path: pathlib.Path, nodes: List[Tuple[float, float, float]], elements: List[List[int]]) -> None:
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with file_path.open("w", encoding="utf-8") as out:
        out.write("<?xml version=\"1.0\"?>\n")
        out.write("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
        out.write("  <UnstructuredGrid>\n")
        out.write(f"    <Piece NumberOfPoints=\"{len(nodes)}\" NumberOfCells=\"{len(elements)}\">\n")
        out.write("      <PointData Scalars=\"NodeId\">\n")
        out.write("        <DataArray type=\"Int32\" Name=\"NodeId\" format=\"ascii\">\n")
        for i in range(len(nodes)):
            out.write(f"          {i + 1}\n")
        out.write("        </DataArray>\n")
        out.write("      </PointData>\n")
        out.write("      <Points>\n")
        out.write("        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n")
        for x, y, z in nodes:
            out.write(f"          {x} {y} {z}\n")
        out.write("        </DataArray>\n")
        out.write("      </Points>\n")
        out.write("      <Cells>\n")
        out.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
        for element in elements:
            out.write("          " + " ".join(str(node) for node in element) + "\n")
        out.write("        </DataArray>\n")
        out.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
        offset = 0
        for _ in elements:
            offset += 6
            out.write(f"          {offset}\n")
        out.write("        </DataArray>\n")
        out.write("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")
        for _ in elements:
            out.write("          22\n")
        out.write("        </DataArray>\n")
        out.write("      </Cells>\n")
        out.write("    </Piece>\n")
        out.write("  </UnstructuredGrid>\n")
        out.write("</VTKFile>\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert a quadratic Gmsh .msh mesh to .vtu for ParaView.")
    parser.add_argument("input_msh", type=pathlib.Path)
    parser.add_argument("output_vtu", type=pathlib.Path)
    args = parser.parse_args()

    nodes, elements = parse_msh_ascii_v4(args.input_msh)
    write_vtu(args.output_vtu, nodes, elements)


if __name__ == "__main__":
    main()
