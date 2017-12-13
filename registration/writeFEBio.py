# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
from scipy.stats.mstats import find_repeats
import pprocess

def NodeDisp(NC,Elem,nodeConst,nodes,disp,filename):
  fid = open(filename+'.feb','w')
  fid.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="1.0">\n')
  fid.write('	<Control>\n')
  fid.write('		<title>'+filename+'</title>\n')
  fid.write('		<time_steps>5</time_steps>\n')
  fid.write('		<step_size>0.2</step_size>\n')
  fid.write('		<max_refs>15</max_refs>\n')
  fid.write('		<max_ups>10</max_ups>\n')
  fid.write('		<dtol>0.001</dtol>\n')
  fid.write('		<etol>0.01</etol>\n')
  fid.write('		<rtol>0</rtol>\n')
  fid.write('		<lstol>0.9</lstol>\n')	
  fid.write('		<pressure_stiffness>1</pressure_stiffness>\n')
  fid.write('		<time_stepper>\n')
  fid.write('			<dtmin>0.01</dtmin>\n')
  fid.write('			<dtmax>0.1</dtmax>\n')
  fid.write('			<max_retries>10</max_retries>\n')
  fid.write('			<opt_iter>10</opt_iter>\n')
  fid.write('		</time_stepper>\n')
  fid.write('		<plot_level>PLOT_DEFAULT</plot_level>\n')
  fid.write('	</Control>\n')
  fid.write('	<Material>\n')
  fid.write('		<material id="1" name="isoBone" type="linear elastic">\n')
  fid.write('			<E>1.6e+006</E>\n')
  fid.write('			<v>0.3</v>\n')
  fid.write('		</material>\n')
  fid.write('	</Material>\n')
  fid.write('	<Geometry>\n')
  fid.write('		<Nodes>\n')
  for i in range(NC.shape[0]):
    fid.write('			<node id="'+np.str(i+1)+'"> '+np.str(NC[i,0])+', '+np.str(NC[i,1])+', '+np.str(NC[i,2])+'</node>\n')
  fid.write('		</Nodes>\n')
  fid.write('		<Elements>\n')
  Elem = Elem+1
  for i in range(Elem.shape[0]):
    fid.write('			<tet4 id="'+np.str(i+1)+'" mat="1">'+np.str(Elem[i,0])+','+np.str(Elem[i,1])+','+np.str(Elem[i,2])+','+np.str(Elem[i,3])+'</tet4>\n')
  fid.write('		</Elements>\n')
  fid.write('	</Geometry>\n')
  fid.write('	<Boundary>\n')
  if nodeConst.size>0:
    fid.write('		<fix>\n')
    if nodeConst.size>1:
      for i in range(nodeConst.size):
	fid.write('			<node id="'+np.str(nodeConst[i]+1)+'" bc="xyz"></node>\n')
    else:
      fid.write('			<node id="'+np.str(nodeConst+1)+'" bc="xyz"></node>\n')
      fid.write('			<node id="'+np.str(nodeConst+1)+'" bc="uvw"></node>\n')
    fid.write('		</fix>\n')
  fid.write('		<prescribe>\n')
  for i in range(nodes.size):
    fid.write('			<node id="'+np.str(nodes[i]+1)+'" bc="x" >'+np.str(disp[i,0])+'</node>\n')
    fid.write('			<node id="'+np.str(nodes[i]+1)+'" bc="y" >'+np.str(disp[i,1])+'</node>\n')
    fid.write('			<node id="'+np.str(nodes[i]+1)+'" bc="z" >'+np.str(disp[i,2])+'</node>\n')
  fid.write('		</prescribe>\n')
  fid.write('	</Boundary>\n')
  fid.write('	<Output>\n')
  fid.write('		<logfile file="'+filename+'NC.txt">\n')
  fid.write('			<node_data data="x;y;z"></node_data>\n')
  fid.write('		</logfile>\n')
  fid.write('	</Output>\n')
  fid.write('</febio_spec>\n')
  fid.close()
  
  
def NodeForce(NC,Elem,nodeConst,nodeLists,loadC,filename):
  fid = open(filename+'.feb','w')
  fid.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="1.0">\n')
  fid.write('	<Control>\n')
  fid.write('		<title>'+filename+'</title>\n')
  fid.write('		<time_steps>1</time_steps>\n')
  fid.write('		<step_size>1</step_size>\n')
  fid.write('		<max_refs>15</max_refs>\n')
  fid.write('		<max_ups>10</max_ups>\n')
  fid.write('		<dtol>0.001</dtol>\n')
  fid.write('		<etol>0.01</etol>\n')
  fid.write('		<rtol>0</rtol>\n')
  fid.write('		<lstol>0.9</lstol>\n')	
  fid.write('		<pressure_stiffness>1</pressure_stiffness>\n')
  fid.write('		<time_stepper>\n')
  fid.write('			<dtmin>0.01</dtmin>\n')
  fid.write('			<dtmax>0.1</dtmax>\n')
  fid.write('			<max_retries>10</max_retries>\n')
  fid.write('			<opt_iter>10</opt_iter>\n')
  fid.write('		</time_stepper>\n')
  fid.write('		<plot_level>PLOT_DEFAULT</plot_level>\n')
  fid.write('	</Control>\n')
  fid.write('	<Material>\n')
  fid.write('		<material id="1" name="isoBone" type="linear elastic">\n')
  fid.write('			<E>1.6e+006</E>\n')
  fid.write('			<v>0.3</v>\n')
  fid.write('		</material>\n')
  fid.write('	</Material>\n')
  fid.write('	<Geometry>\n')
  fid.write('		<Nodes>\n')
  for i in range(NC.shape[0]):
    fid.write('			<node id="'+np.str(i+1)+'"> '+np.str(NC[i,0])+', '+np.str(NC[i,1])+', '+np.str(NC[i,2])+'</node>\n')
  fid.write('		</Nodes>\n')
  fid.write('		<Elements>\n')
  Elem = Elem+1
  for i in range(Elem.shape[0]):
    fid.write('			<tet4 id="'+np.str(i+1)+'" mat="1">'+np.str(Elem[i,0])+','+np.str(Elem[i,1])+','+np.str(Elem[i,2])+','+np.str(Elem[i,3])+'</tet4>\n')
  fid.write('		</Elements>\n')
  fid.write('	</Geometry>\n')
  fid.write('	<Boundary>\n')
  if nodeConst.size>0:
    fid.write('		<fix>\n')
    if nodeConst.size>1:
      for i in range(nodeConst.size):
	fid.write('			<node id="'+np.str(nodeConst[i]+1)+'" bc="xyz"></node>\n')
    else:
      fid.write('			<node id="'+np.str(nodeConst+1)+'" bc="xyz"></node>\n')
      fid.write('			<node id="'+np.str(nodeConst+1)+'" bc="uvw"></node>\n')
    fid.write('		</fix>\n')
  fid.write('		<force>\n')
  for i in range(len(loadC)):
    for j in nodeLists[i]:
      fid.write('			<node id="'+np.str(j+1)+'" bc="'+loadC[i][0]+'" lc="'+np.str(loadC[i][1])+'">'+np.str(loadC[i][2])+'</node>\n')
  fid.write('		</force>\n')
  fid.write('	</Boundary>\n')
  fid.write('	<LoadData>\n')
  for i in range(len(loadC)):
    fid.write('		<loadcurve id="'+str(loadC[i][1])+'">\n')
    fid.write('			<loadpoint>0,0</loadpoint>\n')
    fid.write('			<loadpoint>1,'+str(loadC[i][3])+'</loadpoint>\n')
    fid.write('		</loadcurve>\n')
  fid.write('	</LoadData>\n')
  fid.write('	<Output>\n')
  fid.write('		<logfile file="'+filename+'NC.txt">\n')
  fid.write('			<element_data data="sx;sy;sz" name ="element stress x, y and z"></element_data>\n')
  fid.write('			<element_data data="sxy;syz;sxz" name ="element stress xy, yz and xz"></element_data>\n')
  fid.write('		</logfile>\n')
  fid.write('	</Output>\n')
  fid.write('</febio_spec>\n')
  fid.close()