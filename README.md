# SPP2231
Implementations used for my PhD time in the FluSimPro SPP 2231 project. The work is mostly related to spray cooling and is currently implemented in OpenFOAM 12.

# Implementations and functionalities
- General injection model "GeneralCoupledSizeVelocityInjection" which can inject parcels with correlated sizes and velocities according to normal distributions or cumulative distribution functions (CDF). More details can be found in: 

https://www.researchgate.net/publication/395319184_EulerLagrange_simulation_of_different_injection_models_based_on_PDA_measurements 

Standard OpenFOAM now has something similar for the distribution and such a structure could be used (this version was adapted from my older version in OF8 and was not adapted).

- Dispersion model "StochasticLangevinDispersionRAS" which can be found in the Best Practices Guidelines from Prof. Martin Sommerfeld. More details and description will be added.

- Relevant time steps based on particle's response time and turbulent dispersion.

- Particle collector to register specific information of parcels crossing defined planes.

# Authors
- Kaissar de Oliveira Nabbout

https://www.researchgate.net/profile/Kaissar-Nabbout

https://www.linkedin.com/in/kaissar-nabbout/

# License
This project is licensed under the GPL v3 License - see the LICENSE.md file for details.
