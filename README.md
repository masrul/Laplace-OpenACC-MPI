# Hybrid GPU/CPU programing for Laplacain 

This repo contains OpenACC/MPI version of laplacian solution. This was a hackathon problem in International HPC summer schol at Slovenia. I was one of the contestants. My solution postioned 2nd best timing in hybrid GPU/CPU categories. Solutions was tested in   Bridge machine of Pittsburgh Supercomputing Center (PSC). I used 4 computing nodes (maximum allowed). Each node has 4 GPUs. So i launched 16 MPI processes, each for each GPU. 

My solution took ``7.53 sec``, which  scored second postion. The best timing was ``~7.35``  sec. 
