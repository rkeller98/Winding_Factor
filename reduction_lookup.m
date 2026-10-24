function y = reduction_lookup(map,x)
z=reduction_coordinates(x,map.model,map.p,map.current_scale);
z=(z-map.origin)./map.scale; y=nan(size(x,1),2);
args=num2cell(z,1);
for c=1:2
    F=griddedInterpolant(map.axes,map.values{c},'linear','none');
    y(:,c)=F(args{:});
end
end
