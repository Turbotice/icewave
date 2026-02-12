function d = bound_harmonicN_fit(x,y,N,h_w)
    yth = bound_harmonicN(x,N,h_w);
    d = sum((yth - y).^2);
end