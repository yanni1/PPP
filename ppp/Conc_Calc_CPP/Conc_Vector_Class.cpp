#include <vector>

class ConcArray {
    private:
        int species_count, ny, nx;
        std::vector<float> data;
    
    public:
        ConcArray(int species_count, int ny, int nx)
            : species_count(species_count), ny(ny), nx(nx),
              data(species_count * ny * nx, 0.0f) {}
    
        // Flat index calculation
        inline int index(int s, int y, int x) const {
            return (s * ny + y) * nx + x;
        }
    
        // Accessor
        float& operator()(int s, int y, int x) {
            return data[index(s, y, x)];
        }
    
        const float& operator()(int s, int y, int x) const {
            return data[index(s, y, x)];
        }
    
        // Getter for raw vector (e.g., if needed by external functions)
        std::vector<float>& raw() { return data; }
        const std::vector<float>& raw() const { return data; }
    
        // Resize or reset methods could be added as needed
    };
    