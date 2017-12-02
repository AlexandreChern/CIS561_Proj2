using FactCheck
using PETSc2


facts("  ---Checking Petsc data types---") do

  @fact PetscScalar => Float64
  @fact PetscReal => Float64
  @fact PetscInt => Int32

end

FactCheck.exitstatus()
