using Base.Test
using PETSc2

@testset "  ---Checking Petsc data types---" begin

  @test ( PetscScalar )== Complex128
  @test ( PetscReal )== Float64
  @test ( PetscInt )== Int32

end


