function H=heatkernel(V,v, t)
  H = V * exp(t*diag(v)) * V';
end
