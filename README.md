# tadpole-spinal-cord-models

This is the source code for three compuational models used for generating connectivity and swimming functionality of spinal cord neurons in the Xenopus 
tadpoles using biological data. The model is an update of previous models models of the Xenopus tadpole swimming spinal network (Roberts et al. 2014; 
Ferrario et al, 2018) to incorporate new neural populations that initiate and stop locomotor activity. 

1. The first model uses biological data to reconstruct a complete map of synaptic interactions between spinal neuron (anatomical 
connectome), in the folder "anatomical model". 
2. The second model generalizes the anatomical model by averaging a number of anatomical connectomes and thus generating a matrix 
of connection probabilities between the neurons, in the folder "probabilistic model".
3. The third model uses the previous two to generate connectivity between neurons and simulates the dynamics of this circuit 
using Neuron (tested with version 7.3), in the folder "functional model".

A README file contained in the folder of each piece of the model explains how to run the simulations and give a description of the models. 
For a more complete description of these models have a look at (please cite, if you are using this model):

Ferrario, Andrea, et al. "Whole animal modelling reveals neuronal mechanisms of decision-making and reproduces unpredictable swimming in frog tadpoles." bioRxiv (2021).

For any problem contact me. 

Author: Andrea Ferrario, University of Exeter
Last Update: 29/09/2021
Email: A.A.Ferrario@exeter.ac.uk






