@echo off


rem If the R_SNOW_RSCRIPT_CMD variable is not defined assume Rscript
rem is in the search path.

if not defined R_SNOW_RSCRIPT_CMD set R_SNOW_RSCRIPT_CMD=Rscript


%R_SNOW_RSCRIPT_CMD% %R_SNOW_LIB%\snow\%*
