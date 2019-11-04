from ij import IJ, WindowManager, Prefs
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import Roi, Toolbar, WaitForUserDialog
import os
from ij.io import DirectoryChooser, OpenDialog
from ij.process import ImageProcessor
from ij.plugin import ChannelSplitter

# Bioformats
from loci.plugins import BF
from loci.common import Region

def close_images():
	open_images = WindowManager.getImageTitles()
	for imagename in open_images:
		print imagename
		IJ.selectWindow(imagename)
		IJ.run("Close")
	print "Windows closed"

def nuclei_processor(imp_particles, thresh_type, folder, impname, channel_name):
	rm = RoiManager.getInstance()
	if not rm:
  		rm = RoiManager()
	rm.reset()

	# define a results table
	rt = ResultsTable.getResultsTable()

	imp_particles.show()
	# generate thresholded image for ROI
	nuclei = imp_particles.duplicate()
	nuclei.show()
	IJ.run("Gaussian Blur...", "sigma=3")
	IJ.setAutoThreshold(nuclei, thresh_type)
	IJ.run("Convert to Mask")
	IJ.run("Fill Holes")
#	IJ.run("Watershed")

	# select thresholded image (to set ROIs)
	IJ.run("Set Measurements...", "area mean standard min area_fraction limit display add redirect=["+imp_particles.title+"] decimal=3")
	IJ.run("Analyze Particles...", "size=30-Infinity show=Outlines display clear add")

	# get the ROI manager and save
	rm.runCommand("save selected", os.path.join(folder, impname+'_'+channel_name+"_ROIs.zip"))

#	pixel_collector(rm, imp_measure, channel_name, impname, folder)
#	return rm

def cell_processor(imp_particles, imp_measure, folder, impname, channel_name):

	# Get ROI manager instance, save to zip
	rm = RoiManager.getInstance()
	if not rm:
	  rm = RoiManager()
	rm.reset()

	imp_particles.show()

	# Wait for user to add ROIs to manager
	Prefs.multiPointMode = True
	pause = WaitForUserDialog("Select ROIs of interest and add to manager \n \nPress OK to continue")
	pause.show()

	# get the ROI manager and save
	rm = RoiManager.getInstance()
	rm.runCommand("save selected", os.path.join(folder, impname+"_cells_ROIs.zip"))

	# select image to measure and show
	imp_measure.show()
	pixel_collector(rm, imp_measure, channel_name, impname, folder)

	return rm


def pixel_collector(rm, channel_imp, channel_name, impname, folder):

		# define new Results table
		rt = ResultsTable()

		IndRois = rm.getIndexes()
		for index in IndRois:
			ROI = rm.getRoi(index)
			ROI_name = ROI.getName()
			coords = ROI.getContainedPoints()

			row = 0
			for pixel in coords:
				x_coord = pixel.getX()
				y_coord = pixel.getY()

				rt.setValue(ROI_name+"_X_pos", row, int(x_coord))
				rt.setValue(ROI_name+"_Y_pos", row, int(y_coord))

				pixel_2 = channel_imp.getProcessor().getPixel(int(x_coord), int(y_coord))
				rt.setValue(ROI_name+"_"+channel_name, row, pixel_2)

		  		row = row + 1
		rt.show("Results")

		rt.save(os.path.join(folder, impname+'_'+channel_name+"_pixels.csv"))
		print "Pixel collection done!"

## ------------- processing happens here -----------------
# get input path for merge file
opener = DirectoryChooser("Select the input folder")
input_folder = opener.getDirectory()

#Define results output folder
folder = input_folder+"Results\\"
if not os.path.exists(folder):
	os.mkdir(folder)
# collect list of files to be processed
file_list = [filename for filename in os.listdir(input_folder) if ".tif" in filename]

print file_list

# Define ROI manager, clear if already exists
rm = RoiManager.getInstance()
if not rm:
		rm = RoiManager()
rm.reset()


# Process images
for impfile in file_list:
	# open image
	imp = IJ.openImage(os.path.join(input_folder, impfile))
	imp.show()

	impname=impfile.strip('.tif')
	imp_details = impname.split('_')
	stain, mutant, pos_number = imp_details

	print mutant, stain, pos_number

	imps = ChannelSplitter.split(imp)
	
	# Process cell ROIs
	rm = cell_processor(imp_particles=imps[3], imp_measure=imps[2], folder=folder, impname=impname, channel_name='GFP_cells')
	pixel_collector(rm, imps[3], 'cyto_mCherry_cells', impname, folder)
	pixel_collector(rm, imps[4], 'inc_mCherry_cells', impname, folder)
	rm.reset()

	## Process nuclei image
	imps_thresh = ChannelSplitter.split(imp)
	nuclei_processor(imp_particles=imps_thresh[0], thresh_type='Otsu dark', folder=folder, impname=impname, channel_name='nuclei')
	pixel_collector(rm, imps[2], 'GFP_nuclei', impname, folder)
	pixel_collector(rm, imps[3], 'cyto_mCherry_nuclei', impname, folder)
	pixel_collector(rm, imps[4], 'inc_mCherry_nuclei', impname, folder)
	
	close_images()


print "Process Complete"
print "Results for "+mutant+'_'+pos_number+". Saved in directory: "+folder
