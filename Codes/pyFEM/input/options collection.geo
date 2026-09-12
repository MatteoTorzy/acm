%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
how to cretae deformed views
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fmid,'%s\n','View[0].Visible=0;');
fprintf(fmid,'%s\n','Plugin(Scal2Vec).ViewX=0;');
fprintf(fmid,'%s\n','Plugin(Scal2Vec).ViewY=1;');
fprintf(fmid,'%s\n','Plugin(Scal2Vec).Run;');
fprintf(fmid,'%s\n','View[5].Name = "Def";');
fprintf(fmid,'%s\n','View[5].Visible = 1;');
fprintf(fmid,'%s\n','View[5].IntervalsType=3;');
fprintf(fmid,'%s\n','View[5].VectorType = 5;');
fprintf(fmid,'%s\n','View[5].NbIso=21;');

View[0].Visible=0;
Plugin(Scal2Vec).ViewX=0;
Plugin(Scal2Vec).ViewY=1;
Plugin(Scal2Vec).Run;
View[5].Name = "Def";
View[5].Visible = 1;
View[5].IntervalsType=3;
View[5].VectorType = 5;
View[5].NbIso=21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
how to extract single views from vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mesh.SurfaceEdges=0;

View[0].VectorType = 5;  // (1=segment, 2=arrow, 3=pyramid, 4=3D arrow, 5=displacement, 6=comet)
View[0].Visible = 0;   // Ux 

Plugin(MathEval).View=0;
Plugin(MathEval).Expression0="v0";
Plugin(MathEval).Run;
View[2].Name = "U1";
View[2].Visible = 1;
View[2].NbIso=20;
View[2].IntervalsType=3;

Plugin(MathEval).View=0;
Plugin(MathEval).Expression0="v1";
Plugin(MathEval).Run;
View[3].Name = "U2";
View[3].Visible = 0;
View[3].NbIso=20;
View[3].IntervalsType=3;

View[1].VectorType = 5;  // (1=segment, 2=arrow, 3=pyramid, 4=3D arrow, 5=displacement, 6=comet)
View[1].Visible = 0;   // stress 

Plugin(MathEval).View=1;
Plugin(MathEval).Expression0="v0";
Plugin(MathEval).Run;
View[4].Name = "S11";
View[4].Visible = 0;
View[4].NbIso=20;
View[4].IntervalsType=3;

Plugin(MathEval).View=1;
Plugin(MathEval).Expression0="v1";
Plugin(MathEval).Run;
View[5].Name = "S22";
View[5].Visible = 0;
View[5].NbIso=20;
View[5].IntervalsType=3;

Plugin(MathEval).View=1;
Plugin(MathEval).Expression0="v2";
Plugin(MathEval).Run;
View[6].Name = "S12";
View[6].Visible = 0;
View[6].NbIso=20;
View[6].IntervalsType=3;

