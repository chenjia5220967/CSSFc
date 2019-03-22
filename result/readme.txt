notice: 
This version is for reference ONLY as it is reorganized from the raw version without test.
For better performance of this code, the sequence images should be as similar as possible.

file structure:
	train: code for training
	predict: code for prediction
	common: auxiliary code for training and prediction
	data: MODIS and LandSat-7 images
	model: result of training
	result: fused images

Before running the program, you should download and compile the SPAMS software, and add the "build" folder into PATH.

train script: train/dict_train
prediction script: predict/CSDL_predict