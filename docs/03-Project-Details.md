# **3** Project Details {-#project-details}

## **3.1** Objective {-#objective}

The objective of this project is to develop open source tools to integrate current ALS-based forest inventory attributes with pre-existing Ontario-wide Forest Resources Inventory (FRI) data to map key forest attributes of interest (specifically age and species composition) in a polygonal format over entire forest management areas in Ontario. 

## **3.2** Themes {-#themes}

* Enhance the production of the Enhanced Forest Resource Inventory (eFRI)

* Integrate and refine photo interpreted attributes and ALS

* Develop object-based image assessment techniques to segment forest stands from ALS attributes

* Build a methodology to impute age and species composition from the FRI into newly-segmented forest stands

* Deliver open source code and methods

## **3.3** Description {-#description}

A key requirement of sustainable forest management is the establishment and maintenance of forest inventories to provide accurate and timely information on the state of the forest to support a variety of purposes and information needs. In Ontario, where extensive forest management practices dominate, forest inventories are derived from a multi-stage process that involves acquiring aerial photography, using the photos to delineate homogenous units or forest stands, and then interpreting attributes from the photography for those delineated stands, often with the aid of stereo vision. 

Due to the capability of ALS to provide a detailed three-dimensional view of forest structure, the technology and related analysis approaches have transformed the derivation of forest inventories. Attributes such as volume, biomass, diameter, and basal area as well as information on canopy cover and its vertical distribution have all been shown to meet or exceed accuracy requirements in operational forest management with relative root mean square errors (%RMSE) typically lower than 20%.  

A key issue in utilising ALS data to derive forest inventory information is that a number of forest attributes required to provide the full complement of variables necessary for strategic-level forest inventories such as species group, composition, stand age, and site productivity, are not well predicted from ALS data alone, due to the lack of spectral and other information in the returned point cloud. A second issue is that the development of Enhanced Forest Inventories (EFI) using the area-based approaches (ABA) results in raster-based predictions of forest attributes. While these raters are produced at a fine grid cell (typically 20 x 20 m), most forest management decision making processes, such as harvest planning and scheduling, growth and yield estimation and longer term estate planning relies on polygonal data structures, such as those manually derived from photo interpreters.

As a result, two gaps are evident when introducing ALS-based inventories into a forest management framework for strategic and tactical decision making. First, raster based EFI information needs to be transformed into polygonal structures and second, new approaches are needed to derive the attributes required for sustainable forest management decision making that are not readily predicted from ALS data alone.

The goal of this project is to develop and apply an open source methodological framework that combines raster ALS/EFI predictions of stand structure with environmental and satellite data to produce comprehensive polygonal forest inventory information in a spatially exhaustive way, over large forested areas of Ontario. These new layers will be polygon based, with the full set of required attributes, at a scale comparable to existing photo-interpreted inventory, to ensure they are compatible with subsequent forest growth and yield, harvest and estate planning software and activities. Once new forest stand polygons are generated, conventional FRI attributes such as species and age are carried forward and integrated using an imputation algorithm.

## **3.4** Study Area {-#study-area}

We propose two forest management units as focus sites for this research.

The Romeo Malette Forest, which is an actively managed forest by GreenFirst Forest Products and is a typical example of a boreal forest management unit in this region. Forest harvesting practices result in an array of stand development stages which are associated with a range of vegetation structures. 

A second focus site will be French Severn forest. This area has a distinct and complementary ecology, resulting in a range of species and structures not evident in the northern boreal.

## **3.5** Methodology {-#methodology}

First, we employ image segmentation (geographic object-based image analysis) to generate new up-to-date forest stands using key ALS-derived forest attributes of forest height, canopy cover, and coefficient of variation.

Second, we ask what is the optimal selection of primary ALS and modelled EFI attributes to drive the imputation of photo-interpreted forest composition attributes, and what is the impact of imputation parameters on the correspondence between imputed and interpreted forest age and species composition. 

After establishing an optimal imputation model, we then undertake imputation to examine how the distribution of key forest attributes changes from the conventional FRI to an integrated inventory containing both ALS and imputed photo-interpreted attributes.

Methodological details along with code demonstrations and results are contained in the appropriate sections below.

This work is guided by the need to present the forest management community with an open-source, reproducible workflow to update forest inventories, while maintaining the forest stand as the baseline unit and integrating ALS data with important attributes still contained in conventional photo-interpreted forest inventories.

## **3.6** Knowledge & Technology Transfer {-#knowledge-transfer}

This project is designed to support research, development and technology transfer in the use of transformation technologies such as innovative remote sensing and environmental datasets such as ALS for forest management across Canada’s forest sector. Specifically this project is designed to support the Provincial government of Ontario to improve the accuracy of the forest inventory through the innovative use of ALS data, as well as exploiting the available information on species and age in the current FRI. By undertaking these activities we contribute to the ongoing transformation of the forest sector through the development and adoption of innovative science-based solutions in particular by linking to the new Ontario governments Forest Sector Strategy. 

As the application of ALS and other 3-dimensional data to forest inventories becomes more common, the need for joint projects with industry is required to communicate the successes and caveats of these technologies and their appropriate application across Canada. The ultimate outcomes of this project will be a workforce that is more informed on the use of remote sensing technologies for next generation forest inventory applications.

We will work with Ontario FRI staff and GreenFirst Forest Products managers on a number of key methods of outcome dissemination. These will include:

* Peer-reviewed scientific publications

* A series of online workshops

* A technical report describing the methods (this document)

* Open access to R and other developed programming scripts to undertake the modelling and comparisons. 

More broadly key outputs of this project will be an increased awareness and capacity building of operational forest managers on the use of these technologies for next generation forest inventory which will lead to an increased ability to dissolve current boundaries between operational (tactical) and strategic (forest level) inventories.

In addition, more thorough adoption of ALS technologies as a critical component of forest inventories will benefit other Canadian industries such as survey providers and technology developers by increasing the market and consolidating their international competitiveness, key goals of the Ontario Forest Sector strategy.
