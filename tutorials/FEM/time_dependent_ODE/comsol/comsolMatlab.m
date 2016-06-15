%function out = model
%
% comsolMatlab.m
%
% Model exported on Jun 2 2016, 23:04 by COMSOL 5.2.0.166.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/giacomo/Google Drive/MODEL/tutorials/FEM/time_dependent_ODE/comsol');

model.comments(['untitled\n\n']);

model.modelNode.create('comp1');

model.geom.create('geom1', 3);

model.mesh.create('mesh1', 'geom1');

model.physics.create('c', 'CoefficientFormPDE', 'geom1', {'u'});

model.geom('geom1').create('blk1', 'Block');
model.geom('geom1').feature('blk1').set('base', 'center');
model.geom('geom1').run('blk1');
model.geom('geom1').run;

model.physics('c').feature('cfeq1').setIndex('c', {'0' '0' '0' '0' '0' '0' '0' '0' '0'}, 0);
model.physics('c').feature('cfeq1').setIndex('a', '1', 0);
model.physics('c').feature('cfeq1').setIndex('f', '(x==0)*(y==0)*(z==0)', 0);
model.physics('c').feature('cfeq1').setIndex('da', '0', 0);

model.mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study('std1').feature('stat').activate('c', true);
model.study('std1').feature('stat').set('showdistribute', false);

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');

model.result.create('pg1', 3);
model.result('pg1').set('data', 'dset1');
model.result('pg1').create('slc1', 'Slice');
model.result('pg1').feature('slc1').set('expr', 'u');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').run;

model.physics('c').feature('cfeq1').setIndex('f', 'x^2+y^2+z^2<0.01', 0);

model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').feature.remove('s1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;

model.physics('c').feature('cfeq1').setIndex('f', 'x^2+y^2+z^2<0.001', 0);

model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').feature.remove('s1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;

model.study('std1').feature('stat').set('showdistribute', false);

model.physics('c').prop('ShapeProperty').set('order', '1');

model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').feature.remove('s1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;

model.physics('c').feature.create('weak1', 'WeakContribution', 3);
model.physics('c').feature('weak1').selection.set([1]);
model.physics('c').feature.create('weak2', 'WeakContribution', 0);
model.physics('c').feature('weak2').selection.set([6]);
model.physics('c').feature('weak2').set('weakExpression', '1');
model.physics('c').feature.remove('weak1');
model.physics('c').feature('cfeq1').setIndex('f', '0', 0);

model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').feature.remove('s1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;

model.geom('geom1').run('blk1');
model.geom('geom1').create('pt1', 'Point');
model.geom('geom1').run('pt1');
model.geom('geom1').run;

model.physics('c').feature('weak2').selection.set([]);

model.view('view1').set('transparency', 'on');

model.physics('c').feature('weak2').selection.set([5]);

model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').feature.remove('s1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;

model.physics('c').feature('weak2').set('weakExpression', 'test(u)');

model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').feature.remove('s1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;

model.mesh('mesh1').autoMeshSize(8);
model.mesh('mesh1').run;
model.mesh('mesh1').autoMeshSize(9);
model.mesh('mesh1').run;

model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').feature.remove('s1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').set('allowtableupdate', false);
model.result('pg1').set('title', 'Slice: Dependent variable u (1)');
model.result('pg1').set('hasbeenplotted', true);
model.result('pg1').feature('slc1').set('rangeunit', '1');
model.result('pg1').feature('slc1').set('rangecolormin', -4.726449482019564);
model.result('pg1').feature('slc1').set('rangecolormax', 26.32537964366489);
model.result('pg1').feature('slc1').set('rangecoloractive', 'off');
model.result('pg1').feature('slc1').set('rangedatamin', -4.726449482019564);
model.result('pg1').feature('slc1').set('rangedatamax', 26.32537964366489);
model.result('pg1').feature('slc1').set('rangedataactive', 'off');
model.result('pg1').feature('slc1').set('rangeactualminmax', [-4.726449482019564 26.32537964366489]);
model.result('pg1').feature('slc1').set('hasbeenplotted', true);
model.result('pg1').set('renderdatacached', false);
model.result('pg1').set('allowtableupdate', true);
model.result('pg1').set('renderdatacached', true);
model.result.table.create('evl3', 'Table');
model.result.table('evl3').comments('Interactive 3D values');
model.result.table('evl3').label('Evaluation 3D');
model.result.table('evl3').addRow([2.6644797479491444E-12 -0.3121165480320153 -0.23103330803367694 0.21454516544881808], [0 0 0 0]);
model.result.export.create('tbl1', 'Table');
model.result.export('tbl1').set('filename', '/Users/giacomo/Desktop/Untitled.txt');
model.result.export('tbl1').run;
model.result.export.create('data1', 'dset1', 'Data');
model.result.export('data1').set('expr', {'u'});
model.result.export('data1').set('descr', {'Dependent variable u'});
model.result.export('data1').set('unit', {'1'});
model.result.export('data1').set('filename', '/Users/giacomo/Desktop/Untitled.txt');
model.result.export('data1').run;
model.result.export('data1').set('smooth', 'none');
model.result.export('data1').set('filename', '/Users/giacomo/Google Drive/MODEL/tutorials/FEM/time_dependent_ODE/sol.txt');
model.result.export('data1').run;

model.mesh('mesh1').autoMeshSize(8);
model.mesh('mesh1').run;

out = model;

MA=mphmatrix(model,'sol1','out',{'K','D','E','L','uscale'})
save('../MA.mat','MA')

%% Save mesh in MODEL format
vtx=model.mesh('mesh1').getVertex;
tri=model.mesh('mesh1').getElem('tet');

file_V = fopen('../N/N_6.txt','w');
nodeformat='%i %e %e %e \n';
for k=1:size(vtx,2)
    fprintf(file_V,nodeformat, [k-1 vtx(:,k)']);
end
fclose(file_V);

file_T = fopen('../T/T_6.txt','w');
triformat='%i %i %i %i %i %i \n';
for k=1:size(tri,2)
    fprintf(file_T,triformat, [k-1 tri(:,k)' 1]);
end
fclose(file_T);

%% 
figure(1)
clf
plot(MA.L/max(abs(MA.L)))
b=MA.L*0;
b(31)=-1;
x=MA.K\b;
hold on
plot(x/max(abs(x)),'r')
