This package of Matlab codes implement semi-analytical propagator matrix methods to compute surface deformation in layered elastic or viscoelastic halfspaces due to either fault or mogi (point inflationary) sources. 

There are three codes in this package:

1. PLate_over_Maxwell codes (Fault_Source_Viscoelastic folder). View run_postseismic_example.m to see how to execute these codes.

This code computes the velocities or displacements at surface points at a single time due to uniform slip
  on a rectangular dislocation in an elastic plate overlying a Maxwell
  viscoelastic layer that lies over a Maxwell viscoelastic halfspace. 

  2. Layered elastic halfspace codes (Fault_Source_Elastic folder). View run_example.m to see how to execute these codes.

  LayeredGreens.m calculates the surface displacements due to a uniform, rectangular dislocation in a layered
  elastic halfspace using propagator matrix methods.

  3. Mogi souce codes (Mogi_Souce folder). View run_example.m to see how to execute these codes.

     MogiLayers.m is the same implementation of propagator matrix solutions as in LayeredGreens.m (above), but for inflationary point sources (the 'Mogi' source).
     
