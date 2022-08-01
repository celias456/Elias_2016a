function VP = vmap(Pbar,mu,a,b,ro)
% Defines the v-map for a model
% x(t) = a1 *mu* E_t1(x(t+1)) + a1 *(1-mu)* E_t2(x(t+1)) + a2 * x(t-1) + b1 * z(t)
% z(t) = ro * z(t-1) + v(t)

VP = a*mu*Pbar(2)*ro + b;

end