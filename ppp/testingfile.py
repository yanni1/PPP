import Conc_Calc_CPP

print("Nanobind module works!")

# Call the parallel_add function from C++
result = Conc_Calc_CPP.parallel_add(3, 4)
print(f"Result of 3 + 4 = {result}")
