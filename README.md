🚀2D FDTD Electromagnetic Simulator 
👋 simple 2D FDTD (Finite-Difference Time-Domain) simulator in C++.
<img width="1559" height="787" alt="image" src="https://github.com/user-attachments/assets/b9451f70-0855-4271-93ef-af03229fc68c" />

Simulates 2D Electromagnetic Waves: Based on the classic Yee algorithm.
Handles Real Materials: Uses the Drude model for metals like gold or copper (no more boring perfect conductors!).
mart Boundaries: Has CPML (think of it as "invisibility cloaks" for the edges of the simulation) to prevent waves from bouncing back.
Multiple Sources: You can have more than one "light source" and even control their phase (in-phase, out-of-phase, etc.).
Looks Pretty: Outputs .vtk files so you can visualize the results in ParaView , which is awesome.
Runs Faster: Uses OpenMP to spread the work across all CPU cores.
  
the results! 

The simulation will generate a bunch of fdtd_fields_*.vtk files.

paraview fdtd_fields_*.vtk
