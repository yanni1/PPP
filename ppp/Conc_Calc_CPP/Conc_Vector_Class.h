#ifndef CONCARRAY_H
#define CONCARRAY_H

#include <vector>
#include <stdexcept>

class ConcArray {
private:
    int species_count, ny, nx;
    std::vector<float> data;

public:
    // Constructor
    ConcArray(int species_count, int ny, int nx)
        : species_count(species_count), ny(ny), nx(nx),
          data(species_count * ny * nx, 0.0f) {}

    // Disable copy constructor and assignment if needed
    // ConcArray(const ConcArray&) = delete;
    // ConcArray& operator=(const ConcArray&) = delete;

    // Compute flat index
    inline int index(int s, int y, int x) const {
        return (s * ny + y) * nx + x;
    }

    // Element access (mutable)
    inline float& operator()(int s, int y, int x) {
        return data[index(s, y, x)];
    }

    // Element access (const)
    inline const float& operator()(int s, int y, int x) const {
        return data[index(s, y, x)];
    }

    // Access raw data
    inline std::vector<float>& raw() { return data; }
    inline const std::vector<float>& raw() const { return data; }

    // Getters for dimensions
    inline int get_species_count() const { return species_count; }
    inline int get_ny() const { return ny; }
    inline int get_nx() const { return nx; }

    // Optional: clear or reset
    void reset(float value = 0.0f) {
        std::fill(data.begin(), data.end(), value);
    }
};

#endif // CONCARRAY_H
