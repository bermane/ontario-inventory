# Project Details {-#project-details}

## Objective {-#objective}

The objective of this project is to develop open source tools to integrate current EFI SPL-based forest inventory attributes with pre-existing Ontario-wide Forest Resources Inventory data (FRI) to map key forest attributes of interest (volume, height, basal area, density, species composition and group, age and micro-site productivity) in a polygonal format, over entire forest management areas in Ontario. 

## Themes {-#themes}

*Enhance the production of the Enhanced Forest Resource Inventory (eFRI)

*Integrate and refine photo interpreted attributes and LiDAR

*Develop object-based image assesment techniques to use with LiDAR for species identification

*Deliver open source code and methods

## Description {-#description}

A key requirement of sustainable forest management is the establishment and maintenance of forest inventories to provide accurate and timely information on the state of the forest to support a variety of purposes and information needs. In Ontario, where extensive forest management practices dominate, forest inventories are derived from a multi-stage process that involves acquiring aerial photography, using the photos to delineate homogenous units or forest stands, and then interpreting attributes from the photography for those delineated stands, often with the aid of stereo vision. 

Due to the capability of airborne Single Photon LiDAR (SPL) to provide a detailed three-dimensional view of forest structure, the technology and related analysis approaches have transformed the derivation of forest inventories. Attributes such as volume, biomass, diameter, and basal area as well as information on canopy cover and its vertical distribution have all been shown to meet or exceed accuracy requirements in operational forest management with relative root mean square errors (%RMSE) typically lower than 20%.  

A key issue in utilising SPL data to derive forest inventory information is that a number of forest attributes required to provide the full complement of variables necessary for strategic-level forest inventories such as species group, composition, stand age, and site productivity, are not well predicted from SPL data alone, due to the lack of spectral and other information in the returned point cloud. A second issue is that the development of Enhanced Forest Inventories (EFI) using the area-based approaches (ABA) results in raster-based predictions of forest attributes. While these raters are produced at a fine grid cell (typically 20 x 20 m), most forest management decision making processes, such as harvest planning and scheduling, growth and yield estimation and longer term estate planning relies on polygonal data structures, such as those manually derived from photo interpreters.

As a result, two gaps are evident when introducing SPL-based inventories into a forest management framework for strategic and tactical decision making. First, raster based EFI information needs to be transformed into polygonal structures and second, new approaches are needed to derive the attributes required for sustainable forest management decision making that are not readily predicted from SPL data alone.


The goal of this project is to develop and apply an open source methodological framework that combines raster EFI predictions of stand structure, with environmental and satellite data to produce comprehensive polygonal forest inventory information in a spatially exhaustive way, over large forested areas of Ontario. These new layers will be polygon based, with the full set of required attributes, at a scale comparable to existing photo-interpreted inventory, to ensure they are compatible with subsequent forest growth and yield, harvest and estate planning software and activities.

## Design {-#design}

We propose three forest management units as focus sites for this research.

The Romeo Malette Forest, which is an actively managed forest by RYAM Forest Management and is a typical example of a boreal forest management unit in this region. Forest harvesting practices result in an array of stand development stages which are associated with a range of vegetation structures. 

A second focus site will be forest management areas within the Great Lakes forest FMU, including French Severn forest or Algonquin Provincial Park. These forest areas have distinct and complementary ecologies, resulting in a range of species and structures not evident in the northern boreal. Likewise additional datasets to utilise in the model development will be different, resulting in different approaches than the methodology developed for boreal systems. 

If time and resources are available we would apply the approach in the second year at a third site to ensure the approaches has been applied to all three regions. This would provide the Ministry with worked examples to aid tech/knowledge transfer at the regional level and testing operational scalability across Ontario.

## Methodology {-#methodology}

First, image segmentation (or object based image analysis) will be performed on key EFI-derived raster attributes of stand (Lorey’s) height and canopy cover to derive structurally-homogeneous micro-stand objects. We will use newly developed open source segmentation tools, specifically developed for forestry LIDAR data to produce wall-to-wall polygon representations of EFI raster layers. 

Second, the boundaries of these micro-stand objects will be used to extract information from other EFI layers (volume, basal area, crown closure) transforming them from rasters into attributed micro-stand polygons. Microstands will then be cleaned and merged based on size constraints and differences within and between segments (such as between vs within variations in volume and basal area). We will utilise extensive information already produced by the Ontario MNR on Ontario Forest Resources Inventory Photo Interpretation Specifications which state key area, perimeter and quality control standards for polygon delineation to meet Ontario forestry needs. These standards will be used to ensure that the derived polygons from the raster EFI’s will closely resemble the current FRI polygons in terms of area, shape and perimeter distributions. 

In the third step we will develop an imputation-based approach to derive the remaining attributes not readily predicted from LIDAR including age, species group / species composition or forest unit, and site index. We will utilise the existing Ontario FRI polygons within the forest management area. For each FRI polygon we will extract positional (latitude / longitude), terrain, and climatic information from auxiliary 30 m environmental data. We will also access considerable work completed and underway mapping forest age from historical Landsat data and spatial databases of fire, harvesting, and road disturbances compiled as part of the Boreal Disturbance Database. If available, this type of information will also be compiled and attributed into each delineated polygon. We will then develop an imputation reference dataset which links the four desired attributes to both the structural attributes within the FRI (volume, density, basal area and height), as well as the location, terrain and climate data. This will produce a reference library of every combination of structural and compositional conditions. 

In the final step, we will impute the additional forest inventory attributes for each SPL derived micro-stand stand including age, species group / composition, and site index combination for that given stand based on its nearest neighbours in attribute space. We will also work with other researchers funded by KTDD investigating the use of SPL data for additional attribute such as soils prediction to ensure consistency in the use of environmental attributes. 

We will assess accuracy in a number of ways. For delineation of the micro-stands themselves, we will use GIS spatial overlap algorithms to assess degree of spatial coherence between manually-delineated polygons vs the automatically delineated stands micro-stands from the SPL. For imputation, we will assess agreement using independent validation samples of a hold back of FRI polygons and compare the imputation results with the manually interpreted species, age and site index assessments. 

## Schedule {-#schedule}

### Project Dates {-#project-dates}

September 2021 – March 2022: Compile all available SPL data and plot data. Build necessary EFI at each management area. Apply segmentation approaches on rater predictions
MILESTONE: Compiled SPL and EFI datasets over study areas.
March 2022 – August 2022: Refinement of polygon size and shape based on Provincial input. Development of models for imputation of non-structure variables like species. MILESTONE: Segmentation / object based Algorithm refined, Species / composition databases developed.
August 2022 – March 2023: Full model development of non-structure attributes. Imputation approach applied over all sites. Field program undertaken at key sites to ensure accuracy of predictions, visits to unusual stand conditions to verify predictions. 
MILESTONE: Full algorithm testing, Accuracy assessment.
March 2023 – August 2023: Workshops, code demos, open source code packaging / delivery to all partners. Final validation, Final inventory polygon coverages provided to partners. 
MILESTONE: Open source code, workshops, peer reviewed papers.

### Deliverables {-#deliverables}

Key deliverables for this project are:

1. 15 March 2022: Digital Layers Compiled SPL and EFI datasets over study areas. Provided to ministry and industry staff.
2. 15 August 2022: Object based segmentation approach with validation. Code available for broader scale testing and applications.
3. 15 March 2023: Imputation draft paper developed. Imputation code developed ready for testing.
4. 15 August 2023: Open source code release. Workshop for industry and government participants. Peer reviewed papers on approach.

## Knowledge & Technology Transfer {-#knowledge-transfer}

This project is designed to support research, development and technology transfer in the use of transformation technologies such as innovative remote sensing and environmental datasets such as SPL for forest management across Canada’s forest sector. Specifically this project is designed to support the Provincial government of Ontario to improve the accuracy of the forest inventory through the innovative use of SPL data, as well as exploiting the available information on species and structure in the current FRI. By undertaking these activities we contribute to the ongoing transformation of the forest sector through the development and adoption of innovative science-based solutions in particular by linking to the new Ontario governments Forest Sector Strategy. 

As the application of SPL and other 3-dimensional data to forest inventories becomes more common, the need for joint projects with industry is required to communicate these successes and caveats of these technologies and their appropriate application across Canada. The ultimate outcomes of this project will be a workforce that is more informed on the use of remote sensing technologies for next generation forest inventory applications.

We will work with Ontario MNR FRI staff, and RYAM Forest Management managers on a number of key methods of outcome dissemination. These will include:

*Peer-reviewed scientific publications.

*A series of online workshop on the use of methods for generating the micro-stands and the additional attribution approaches. These will be regional, as needed, to ensure the developed approaches are consistent with regionally specific datasets and needs. Workshops will also be designed to focus on industry forest practitioners, as needed.

*A technical memo describing the methods developed for key technical staff as it is these staff that are applying these tools themselves as well as for consulting firms which can also be hired by smaller SFLs to implement solutions.   

*Open access to R and other developed programming scripts to undertake the modeling and comparisons. 

More broadly key outputs of this project will be an increased awareness and capacity building of operational forest managers on the use of these technologies for next generation forest inventory which will lead to an increased ability to dissolve current boundaries between operational (tactical) and strategic (forest level) inventories. 

In addition, more thorough adoption of SPL technologies as a critical component of forest inventories will benefit other Canadian industries such as survey providers and technology developers by increasing the market and consolidating their international competitiveness, key goals of the Ontario Forest Sector strategy.
