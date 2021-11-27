#include <vector>

/**
 * @brief Scalar field on evenly spaced 2D grid
 *
 */
class ScalarField {
   private:
    const std::size_t width;
    const std::size_t height;
    std::vector<double> val;

   public:
    ScalarField();
};
