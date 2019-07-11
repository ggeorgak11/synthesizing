## Synthesizing Training Data for Object Detection in Indoor Scenes.
G. Georgakis, A. Mousavian, A. C. Berg, J. Kosecka, RSS 2017

This is MATLAB code used to generate the synthetic sets for the paper. Set-up to work for kitchen scenes and generates both a blended and a simple superimposed dataset.

[Project Page](https://cs.gmu.edu/~robot/synthesizing.html)

External dependencies
----------------------
1) Modified Poisson (can be found in this link: https://www.mathworks.com/matlabcentral/fileexchange/39438-fast-seamless-image-cloning-by-modified-poisson-equation)
2) toolbox_nyu_depth_v2. Can be found in this link: http://cs.nyu.edu/~silberman/datasets/nyu_depth_v2.html

Optional Dependencies
--------------------
1) Semantic segmentation approach of "Joint Semantic Segmentation and Depth Estimation with Deep Convolutional Networks" or any other semantic segmentation approach that produces labels for indoor scenes. This is needed if you want to use background scenes other than the ones provided on the project webpage.  

Datasets needed
----------------
1) BigBird

a) Download the objects you want to use from: http://rll.berkeley.edu/bigbird/ .

b) Place all of them in the folder object_data/objects/.

c) Create a txt file that contains a list of the used objects (object_data/our_dataset_objects, we provide an example).

2) Refined masks for all objects from the BigBird dataset. Provided in object_data/our_objects_masks/. A few of the masks for which the graphCut could not produce a good segmentation are not included, so have that in mind when reading through the files.
3) NYU background set. 
	-- Download kitchen raw sets from http://cs.nyu.edu/~silberman/datasets/nyu_depth_v2.html . We have included an example scene in folder background_scenes/kitches/ which is already preprocessed with the nyu toolbox and which includes the semantic segmentations. 
	-- Also the rest of the scenes can be downloaded from the project webpage under 'Data'.
	-- Use utils/nyu_preprocess.m to process each video scene (synchronization, alignment, depth filling).
	-- Use the semantic segmentation approach of "Joint Semantic Segmentation and Depth Estimation with Deep Convolutional Networks" to produce 	  segmentations for each frame, and save the files ending in '_seg.mat'. The files should be save the in the same folder as the color and depth images. If you choose to use your own segmentation approach, then make sure to change the useful labels (params.usf_lbls) in the init_params.m. The useful labels should contain categories such as tables, desks, counters etc.

Running the code
-----------------
The main script is the create_synth_set.m. 
All parameters can be found in the init_params.m. 
The script saves both the annotation in xml format ready to be used for (SSD, Faster R-CNN) and in mat files.
