# Porp-Move

Movement model used for simulating harbour porpoise fine scale movement behaviour. Please 
refer to the scientific publication (Nabe-Nielsen et al. 2013) for detailed documentation.

To run model: 
1. Make sure that the raster-data folder is located in the same directory as the Porp-Move 
  1.0 file. 
2: Open the Porp-Move program using NetLogo >=5.0. 
3: Select area (only data for the Great Belt and homogeneous included). 
4: Press landsc-setup, 
5: Select number of porpoises (usually just one), 
6: Press porps-setup and then go.

Note that the size of the landscape cells are 100 x 100 m in the movement model (see Nabe-
Nielsen et al. 2013) but 400 x 400 m in the population model (Nabe-Nielsen et al. 2014). 
In this version of the movement model each patch covers 13 cells. In the population model 
each patch covered one 400 x 400 m cell.

The monitor "memory-move-porp-0" indicates the relative contribution of the two types of 
behaviour to the following 30-min step: the correlated random walk (CRW) move, controlled 
by the satiation memory (previously called working memory) and the spatial memory move 
(controlled by the reference memory).

Start data is based on a particular satellite-tracked porpoise (pttid 4542). Please refer 
to Nabe-Nielsen et al. (2013) for full documentation of the model.


The model was created as part of the project BRIDGES AS BARRIERS PORPOISE MODELLING: DOES 
THE GREAT BELT BRIDGE HINDER MOVEMENT OF HARBOUR PORPOISES IN THE GREAT BELT funded by 
Fehmern Belt A/S

Copyright (C) 2016, Jacob Nabe-Nielsen <jnn@bios.au.dk>

This program is free softwareyou can redistribute it and/or modify it under the terms 
of the GNU General Public License version 2 and only version 2 as published by the Free 
Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program
if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
MA  02110-1301, USA.


The model was developed and tested using NetLogo version 4.1. Development ended 2011-12-05


Literature

Nabe-Nielsen, J., Tougaard, J., Teilmann, J., Lucke, K. & Forchhammer, M.C. (2013) 
How a simple adaptive foraging strategy can lead to emergent home ranges and increased food 
intake. Oikos, 122, 1307–1316.

Nabe-Nielsen, J., Sibly, R.M., Tougaard, J., Teilmann, J. & Sveegaard, S. (2014) Effects of 
noise and by-catch on a Danish harbour porpoise population. Ecological Modelling, 272, 242–251.

