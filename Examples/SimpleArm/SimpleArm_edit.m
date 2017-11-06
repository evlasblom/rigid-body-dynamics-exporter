function Model = SimpleArm_edit(Model)

bw = 0.02;
Model.rigidbody(Model.dof).visual.bodysizemin = -[bw/2,bw/2,bw/2];
Model.rigidbody(Model.dof).visual.bodysizemax = [bw/2,bw/2,0.2];
Model.visual.bodywidth = bw;

end