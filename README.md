## Pharmacophore dynamics

### Scope of the project
Development of a method to accurately predict protein-protein interaction (PPI) in within a short
time (i.e. a few seconds to maybe minutes). There will be a trade-off between accuracy of the 
prediction and speed. Not only are we looking to predict PPI interaction, but more also 
interaction site (paratope-epitope prediction) and interaction mode (alignment and relative 
position of the proteins to each other). Finally, the predictions can be given a score of likelihood
for binding at this site with this configuration, based on a loss function taking into account 
steric and electronic effects. 

Current state of the art: 
- PPI network prediction: Predicts PPI but not the contact site or alignment of proteins to each
    other. Solely the possible interaction of two proteins is predicted. Usually requires a large 
    dataset of known interactions. 
- Paratope-epitope prediction: Recently a lot of deep learning methods were published. They are able
    to screen a large library of molecules in a short period of time, while indicating possible interactions
    and interaction sites. However, they too require a large high quality dataset for training the models. 
    Therefore, one might argue that they are not universally and easily applicable.
- Protein-protein docking: Aims to figure out in which way the proteins are interacting with each other. 
    These experiments are often hard to set up and require extensive expert knowledge. 
- MD-simulation: Simulating the behavior of two proteins in close proximity can often be used to 
    retrieve more information on the interaction of these proteins. Even though the results obtained
    are quite detailed, these simulations may take a long time.  
  
### Tool MUST HAVEs
- easy application without other data requirements
- easy setup with input of say two PDB files (we can make some restrictions there in terms of content)
- decent performance in terms of speed and accuracy (taking into account the trade-off) to existing methods
- applicable for high-throughput screening
- ranking of solutions
- [optional] user interface for analysis
- [optional] user interface for setup of screening
- [optional] user interface for setup of prediction
- [optional] user interface for further more detailed analysis (i.e. MD)

### Resources 
- brainstorming doc: https://docs.google.com/document/d/1jN3kGP3s32hYwgV7EreyDe6fqRPENOUnEqr1nzrfgfg/edit?usp=sharing
- openMM (Python package for custom MD): http://docs.openmm.org/latest/userguide/index.html
- relevant literature: https://docs.google.com/document/d/1DsM0CMgnx2SXJ9afiPaZR56UdCWc2rR5cjWqxXdKos4/edit?usp=sharing