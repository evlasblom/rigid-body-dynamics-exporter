function Model = Upperbody_edit(Model)

bw = 0.02;
Model.rigidbody(4).visual.bodysizemin = -[bw/2,bw/2,bw/2];
Model.rigidbody(4).visual.bodysizemax = [bw/2,0.1,bw/2];
Model.rigidbody(7).visual.bodysizemin = -[bw/2,0.1,bw/2];
Model.rigidbody(7).visual.bodysizemax = [bw/2,bw/2,bw/2];
Model.visual.bodywidth = bw;

end