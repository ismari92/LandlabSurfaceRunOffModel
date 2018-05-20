# Shape files to ASC

## GQIS procedure
1. Open the shape file in QGIS

        Skripsie\Calibration data\Reworked\Shape folders\2528CC_2m.shp

    This shape file is a set of elevation contour lines in vector format (concatenated vectors from point to point to define the contour line). So this file only has a limited number of elevations, not elevation on a regular grid.

1.  Interpolate a grid from the contour lines.  In QGIS use the menu funtion *Raster/Analysis/Grid interpolation*.  The interpolation algorithm interpolates between the contour lines. It seems that the grid is first created with zero values everywhere, and then some of the pixels in the grid are filled in with the contour-line values.  If the new grid has too many pixels, there are large numbers of zero values between grid lines.  The causes trouble later when interpolation is done.  So there is a fine balance between the grid size and the spacing between the grid lines.  

    The setting that worked best for this data set uses power=4, smoothing=12, radius1=0, radius2=0, maxpoints=0, minpoints=0,angle=0, nodata=0.  Grid size 400 and 400.  The DEM so obtained is not very fine resolution, so a DEM  with finer resolution must be obtained elsewhere.  The GDAL command finally used to do the interpolation is as follows, writing to the file `2528CC_2m_DEM.tif`: 

        gdal_grid -zfield ELEVATION -l 2528CC_2m -a invdist:power=4.0:smothing=12.0:radius1=0.0:radius2=0.0:angle=0.0:max_points=0:min_points=0:nodata=0.0 -outsize 400 400 -of GTiff "Skripsie/Calibration data/Reworked/Shape folders/2528CC_2m.shp" "Skripsie/Calibration data/Reworked/Shape folders/2528CC_2m_DEM.tif"

    This now gives an elevation grid of the whole scene. The catchment area is a much smaller part up in the right top corner.

1.  Clip the raster DEM to cut out only the catchment area using *Raster/Extraction/Clipper* to open the *Clip raster by mask layer* dialog box.  Select the DEM `2528CC_2m_DEM.tif` and the mask layer `A2H062 Catchment`. Enter a name for the output file, e.g., `clipped/clipped2.tif`. Set the raster resolution to be the same as the input, or set a new output file resolution.

        gdalwarp -q -cutline "Skripsie/Calibration data/Reworked/Shape folders/A2H062 Catchment.SHP" -tr 0.000624999999999 0.000624999999999 -of GTiff "Skripsie/Calibration data/Reworked/Shape folders/2528CC_2m_DEM.tif" "Skripsie/Calibration data/Reworked/Shape folders/clipped/clipped2.tif"

1. At this point all pixels outside the masked section are set to zero (0) but the grid is still the same size as the original input grid.  Masking only selects a portion of the data set, making the rest of the data set zero, but without changing the size of the grid.  If the grid is displayed in QGIS, the image is mostly white, except in the masked area where the grid values are shades of gray.  Now the image must be cropped to only contain the masked portion.  Follow the instructions shown here:  

    https://gis.stackexchange.com/questions/240641/how-do-you-change-the-extent-of-a-clipped-raster-layer-in-qgis-to-be-the-mbb

    Step 1: Reclassify the area outside of the image from zero to nodata (-99999) then, Step 2: crop to the rectangular area that removes the excess parts of the grid, leaving only the rectangle (called the bounding box) that contains the non-nodata part.  The cropped image was then named `clipped2cropped.tif`.  

    The masked and cropped image is still in GeoTIFF format. It needs to be converted to ASC.

1.  Install the GDAL windows core msi from `http://download.gisinternals.com/sdk/downloads/release-1600-x64-gdal-2-2-3-mapserver-7-0-7/gdal-202-1600-x64-core.msi`

1. Open a command line window (cmd, conemu or whatever).

1. Change directory (cd) to the location where the shape file is saved.

        cd  "\Skripsie\Calibration data\Reworked\Shape folders\clipped"

2. Add the path to the GDAL executables to the system path. This will allow you to run the GDAL tools in the current directory. The following procedure will make the path available only in the current command line window (if you close the command window the additional path disappears again).  The following assumes that the GDAL tools are installed at `C:\Program Files\GDAL`, adapt to your case if necessary.

        SET PATH=%PATH%;"C:\Program Files\GDAL"


4. Convert the GeoTIFF file to an ASC file using `gdal_translate` (http://www.gdal.org/gdal_translate.html). The output formats are defined in (http://www.gdal.org/formats_list.html). In this case use `AAIGrid` to convert to Arc/Info ASCII Grid.

        Usage: gdal_translate [--help-general] [--long-usage]                                             
            [-ot {Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/                                      
                    CInt16/CInt32/CFloat32/CFloat64}] [-strict]                                          
            [-of format] [-b band] [-mask band] [-expand {gray|rgb|rgba}]                              
            [-outsize xsize[%]|0 ysize[%]|0] [-tr xres yres]                                           
            [-r {nearest,bilinear,cubic,cubicspline,lanczos,average,mode}]                             
            [-unscale] [-scale[_bn] [src_min src_max [dst_min dst_max]]]* [-exponent[_bn] exp_val]*    
            [-srcwin xoff yoff xsize ysize] [-epo] [-eco]                                              
            [-projwin ulx uly lrx lry] [-projwin_srs srs_def]                                          
            [-a_srs srs_def] [-a_ullr ulx uly lrx lry] [-a_nodata value]                               
            [-gcp pixel line easting northing [elevation]]*                                            
            [-mo "META-TAG=VALUE"]* [-q] [-sds]                                                        
            [-co "NAME=VALUE"]* [-stats] [-norat]                                                      
            [-oo NAME=VALUE]*                                                                          
            src_dataset dst_dataset                                                                    
    So in this case use this command to create the ASC file:
    
          gdal_translate -of AAIGrid "clipped2cropped.tif" "clipped2cropped.asc"


6. In the command window type `gdalinfo clipped2cropped.asc` to list the ASC file properties

        > gdalinfo clipped2cropped.asc                                         
        Driver: AAIGrid/Arc/Info ASCII Grid                                    
        Files: clipped2cropped.asc                                             
            clipped2cropped.asc.aux.xml                                     
            clipped2cropped.prj                                             
        Size is 62, 53                                                         
        Coordinate System is:                                                  
        GEOGCS["GCS_WGS_1984",                                                 
            DATUM["WGS_1984",                                                  
                SPHEROID["WGS_84",6378137,298.257223563]],                     
            PRIMEM["Greenwich",0],                                             
            UNIT["Degree",0.017453292519943295]]                               
        Origin = (28.208124999999999,-25.758750000000003)                      
        Pixel Size = (0.000625000000000,-0.000625000000000)                    
        Metadata:                                                              
        AREA_OR_POINT=Area                                                   
        Corner Coordinates:                                                    
        Upper Left  (  28.2081250, -25.7587500) ( 28d12'29.25"E, 25d45'31.50"S)
        Lower Left  (  28.2081250, -25.7918750) ( 28d12'29.25"E, 25d47'30.75"S)
        Upper Right (  28.2468750, -25.7587500) ( 28d14'48.75"E, 25d45'31.50"S)
        Lower Right (  28.2468750, -25.7918750) ( 28d14'48.75"E, 25d47'30.75"S)
        Center      (  28.2275000, -25.7753125) ( 28d13'39.00"E, 25d46'31.13"S)
        Band 1 Block=62x1 Type=Float32, ColorInterp=Gray                       
        Min=1344.009 Max=1543.978                                            
        Minimum=1344.009, Maximum=1543.978, Mean=1426.909, StdDev=47.376     
        NoData Value=-99999                                                  
        Metadata:                                                            
            STATISTICS_MAXIMUM=1543.9783935547                                 
            STATISTICS_MEAN=1426.9091115685                                    
            STATISTICS_MINIMUM=1344.0086669922                                 
            STATISTICS_STDDEV=47.375946004946                                  
    
    Note above that the gdalinfo command also reads `clipped2cropped.asc.aux.xml` and `clipped2cropped.prj` to get additional information about the grid.  The last is quite an important file, with the following contents:

        GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]

    This tells us that the ASC file is in WGS84 coordinates and that the units are in degrees (the corner definitions and the pixel size values). So the pixel size is 0.000625 degrees or 10.908e-6 radians.

    When studying the ASC file with a text editor, it is clear that only the selected part `A2H062 Catchment` has valid values in the grid, but the grid values outside the selected part is marked as `NoData`.


1.  The data calculated are lat/lon/elev and must still be converted to (x,y,z). The lan/lon grid has constant angular pitch in latitude and longitude, but this will convert to unequal pitch in (x,y) because of the location on earth.  Converting from the current WGS84 to a Earth-Centred-Earth-Fixed (ECEF) coordinate system as not undertook here. Instead a calibration was done between known points on the terrain and the same points in the grid.  The pitch in one direction was 68.4 m and in the other direction it was 62.5 m. The grid spacing was then set to the mean value of 65 m in both directions.'
