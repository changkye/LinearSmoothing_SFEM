#include "sfem_io.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace nonlinear_esfem_t6::shared
{

void VtuWriter::write(const std::filesystem::path &file_path,
                      const Mesh &mesh,
                      const Eigen::MatrixXd &displacement,
                      const Eigen::MatrixXd &exact_displacement) const
{
    std::filesystem::create_directories(file_path.parent_path());
    std::ofstream out(file_path);
    if (!out)
    {
        throw std::runtime_error("Failed to open VTK output: " + file_path.string());
    }

    out << std::setprecision(16);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << mesh.nodes.rows()
        << "\" NumberOfCells=\"" << mesh.elements.size() << "\">\n";
    out << "      <PointData Vectors=\"Displacement\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < displacement.rows(); ++i)
    {
        out << "          " << displacement(i, 0) << ' ' << displacement(i, 1) << " 0\n";
    }
    out << "        </DataArray>\n";
    if (exact_displacement.rows() == displacement.rows() && exact_displacement.cols() == displacement.cols())
    {
        out << "        <DataArray type=\"Float64\" Name=\"ExactDisplacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int i = 0; i < exact_displacement.rows(); ++i)
        {
            out << "          " << exact_displacement(i, 0) << ' ' << exact_displacement(i, 1) << " 0\n";
        }
        out << "        </DataArray>\n";
    }
    out << "      </PointData>\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.nodes.rows(); ++i)
    {
        out << "          " << mesh.nodes(i, 0) << ' ' << mesh.nodes(i, 1) << " 0\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto &element : mesh.elements)
    {
        out << "          ";
        for (int a = 0; a < 6; ++a)
        {
            out << element[static_cast<std::size_t>(a)];
            if (a < 5)
            {
                out << ' ';
            }
        }
        out << '\n';
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (std::size_t i = 0; i < mesh.elements.size(); ++i)
    {
        offset += 6;
        out << "          " << offset << '\n';
    }
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < mesh.elements.size(); ++i)
    {
        out << "          22\n";
    }
    out << "        </DataArray>\n";
    out << "      </Cells>\n";
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

} // namespace nonlinear_esfem_t6::shared
