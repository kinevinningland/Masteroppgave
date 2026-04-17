# ReSDDP
## Basics
ReSDDP is a research-prototype model for long-term operational planning of hydropower resources based on stochastic dual dynamic programming (SDDP). The model seeks to find the minimum cost dispatch in power systems with a substantial share of hydropower. Its primary output are the operational strategies for hydropower subsystems. The model is designed to coordinate the treatment of detailed hydropower and fine time-resolution by use of spatial decomposition techniques.     

## Principles
The following main principles apply:
* Code: The model is coded in Julia using the JuMP package for formulating the optimization problems.
* Solvers: A variety of optimization solvers can be used to solve the optimization problems.
* Data: The model supports input data following the format of version 10 from the EMPS model. 

## Disclaimer
The code will be under continuous development during the project period (2024-2027). Breaking changes may occur.

## Acknowledgement
The development of ReSDDP has been funded by the Norwegian Research Council in the project "RES100 – Modeling a 100% Renewable Electricity System", project number 344220. The authors gratefully acknowledge the financial support from the user partners Statkraft, Hafslund ECO, Å Energi, Hydro Energi, and Lyse Produksjon.

## Copyright
Copyright (C) 2024, SINTEF Energy Research

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](./LICENSE) for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.