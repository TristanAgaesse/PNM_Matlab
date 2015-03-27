

nIteration=length(nBreakthroughList);
[a,b]=ind2sub([length(thicknessList),length(randomSeedList)],1:nIteration);
thickness = thicknessList(a);


load('RS-00-10.mat')
nBreakthroughList0010=nBreakthroughList;
nInletLinkList0010=nInletLinkList;
thickness0010=thickness;


nBreakthroughList=[nBreakthroughList0010, nBreakthroughList1120, nBreakthroughList2130];
nInletLinkList=[ nInletLinkList0010, nInletLinkList1120, nInletLinkList2130];
thickness=[ thickness0010, thickness1120, thickness2130];


scatter(thickness,nBreakthroughList)