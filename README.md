# N-Body-Project
PHY480 N-Body project

This project aims to simulate the solar system accurately using the N-Body approach. It includes a functional 4th-order predictor-corrector algorithm implemented in Fortran, located in Code/PCTakeTwo.f90.

It is callable from the command line, and will run the experiment for 1 millions years (simulation time). The algorithm works by: 

- First bootstrapping using second-order oribital mechanics, which only predicts the next timesteps.
- Once a history is filled for each body's motion the predictor corrector can function

- It predicts the next position of each body using a predefined 4th order polynomial in velocity and acceleration
- It then corrects this value at the next time-step according to a different predefined 4th order polynomial
- Each body is moved to their corrected positions and velocities
- The error is recorded between the predicted and corrected values, we can use this to update the timestep size whilst keeping the error inside a certain range
- This allows the predictor-corrector algorithm to simulate days over each time-step, rather than seconds

Please note that this project lacks proper documentation, as I was still learning GitHub. However Code/PCTakeTwo.f90 is commented to explain how the algorithm works.

The primary goal of this project is to replicate the Milankovitch cycles, which are observed in geological records. These cycles are caused by perturbations in orbital distance and interactions with Venus. Relevant figures can be found in the "FiguresSecond" folder. Precise simulation is necessary to reproduce these cycles accurately. The project also includes an energy loss diagram that depicts the overall energy lost by the predictor-corrector algorithm, analogous to simulation accuracy in this context. After simulating for 700,000 years, the algorithm experienced a negligible energy loss of 10 millionths of a percent. Unfortunately, the simulation couldn't be extended due to hardware limitations at the time.
