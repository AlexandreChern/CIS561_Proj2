using FactCheck
using PETSc2
facts("  ---Checking Petsc data types---") do

  @fact PetscScalar => Float32
  @fact PetscReal => Float32
  @fact PetscInt => Int32

end

FactCheck.exitstatus()
