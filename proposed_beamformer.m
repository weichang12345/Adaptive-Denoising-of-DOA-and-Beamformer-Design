function [y w B] = proposed_beamformer(matX, N, theta_s_hat, theta_i_hat, a, d, w, L)

a_s = exp(j*pi*sind(theta_s_hat).*d.');
a_i = exp(j*pi*sind(theta_i_hat).*d.');
w_d = w;  % The desired weight vector

% Find a new beamformer so that its beampattern is as close as possible to
% the desired beampattern (array steering) by method of null steering
for l = 1 : L
    u = (eye(N) - a_i(:,l)*a_i(:,l)'/(a_i(:,l)'*a_i(:,l)))*w_d(:,l);  % Projects the desired signal onto the orthogonal complement of interference angle steering vector
    alpha = 1/(u'*a_s(:,l));                                          % Normalization factor in order to maintain distortionless response
    w(:,l) = alpha*u;                                                 % The projected vector with normalization factor. Although the solution is suboptimal, it still outperforms any other approaches
    y(l) = w(:,l)'*matX(:,l);
    for i = 1 : numel(a(1,:))
        B(i,l) = 20*log10(abs(w(:,l)'*a(:,i)));
    end

end