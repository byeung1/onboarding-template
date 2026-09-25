#pragma once

//brings in std::size_t, unsigned int
#include <cstddef>
//like importing arraylkikst class
#include <vector>

// Starter Grid for the 2D heat-diffusion problem.
//
// The evaluation harness uses operator() to set initial conditions and to read
// results; it never touches your internal storage. Keep this interface,
// everything else is yours.
class Grid {
  private:
    std::size_t rows_;
    std::size_t cols_;

    //using stride to ensure that both the number of rows and columns are divisible by 8
    //it is 8 because the CPU loads memory into RAM in cache lines which are 64 bytes
    //each double is 8 bytes, so 64/8 = 8 doubles
    //if there is no padding, the CPU has to perform misaligned loads
    std::size_t stride_;

    //store the grid as one list for efficiency
    //point is given by i * columns + j, indexed at 0
    std::vector<double> arrayOfHeatTiles_;

  public:
    //fill the grid weith zeros
    Grid(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols), stride_((cols + 7) / 8 * 8), arrayOfHeatTiles_(rows * stride_, 0.0) {}

    double& operator()(std::size_t i, std::size_t j) {
      return arrayOfHeatTiles_[i * stride_ + j];
    };
    double  operator()(std::size_t i, std::size_t j) const {
      return arrayOfHeatTiles_[i * stride_ + j];
    };

    //adding getters so the apply_stencil() can access the dimensions
    //no need for setters because apply_stencil() only needs read permissions
    std::size_t rows() const {
      return rows_;
    }
    std::size_t cols() const {
      return cols_;
    }
    std::size_t stride() const {
      return stride_;
    }

    //pass the raw 1d vector to make the stencil loop better
    std::vector<double>& raw_data() {
      return arrayOfHeatTiles_;
    }
    const std::vector<double>& raw_data() const {
      return arrayOfHeatTiles_;
    }
};  

// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
inline void apply_stencil(const Grid& old_grid, Grid& new_grid) {
  std::size_t rows = old_grid.rows();
  std::size_t cols = old_grid.cols();
  std::size_t stride = old_grid.stride();

  const auto& old_tiles = old_grid.raw_data();
  auto& new_tiles = new_grid.raw_data();

  //copy the top and bottom edges
  for(std::size_t j = 0; j < cols; j++) {
    new_tiles[j] = old_tiles[j];
    new_tiles[(rows - 1) * stride + j] = old_tiles[(rows - 1) * stride + j];
  }

  //copy left and right edges
  //skip corners bcs alr did them
  for(std::size_t i = 1; i < rows - 1; i++) {
    new_tiles[i * stride] = old_tiles[i * stride];
    new_tiles[i * stride + cols - 1] = old_tiles[i * stride + cols - 1];
  }

  //apply the heat-spreading formula
  for(std::size_t i = 1; i < rows - 1; i++) {

    std::size_t rowIdx = i * stride;

    for(std::size_t j = 1; j < cols - 1; j ++) {

      std::size_t center = rowIdx + j;
      std::size_t up = center - cols;
      std::size_t down = center + cols;
      std::size_t left = center - 1;
      std::size_t right = center + 1;

      new_tiles[center] = 0.5 * old_tiles[center] + 0.125 * (old_tiles[up] + old_tiles[down] + old_tiles[left] + old_tiles[right]);
    }
  }
};
