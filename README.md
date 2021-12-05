# Manifold-Model-Mesh-Shape
An implementation of the Manifold Model for Shape Deformation outlined in a paper by Freifeld and Black

# Copyright Information
The Geomproc library was authored by Dr. Oliver Van Kiack and provided under an MIT license  

This model is an implementation of the original research by Oren Freifeld and Michael Black in their paper:  
[Lie Bodies: A 3D Representation of Human Shape](https://link.springer.com/chapter/10.1007/978-3-642-33718-5_1)

The meshes used are from the [TOSCA](http://tosca.cs.technion.ac.il/data/toscahires.txt) and [DFAUST](https://dfaust.is.tuebingen.mpg.de/license.html) datasets for academic purposes (course project) and are copyright their original publishers. Clicking either dataset will bring up the relevant license information for said dataset


## Prerequisites 
The model requires python3 and numpy to run, disturbing the file structure for the datasets will almost certainly cause issues as well.

## Running information
Assuming the prerequisites are met, the model can be run from the command line by running   
```python ManifoldMeshModel.py```   
for a Windows system, substitute python3 for unix systems

By default, the model will load the 11 cat meshes from the TOSCA dataset and provide some analytic information from them. 


## More information 
An explanation of some of the theory behind the model can be found in the video [here](https://www.youtube.com/watch?v=FAFZb96mpfI) and information regarding the implementation, results, and current implementation limitations can be found in [this](https://www.youtube.com/watch?v=5hmKpguMGa0) video 


