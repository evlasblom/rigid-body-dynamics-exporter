function Model = SCARA_edit(Model)

bw = Model.visual.bodywidth; % get default width
Model.rigidbody(4).visual.bodysizemin = -[bw/2,bw/2,0.05];
Model.rigidbody(4).visual.bodysizemax = [bw/2,bw/2,0.20];

end