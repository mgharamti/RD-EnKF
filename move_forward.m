function x_new = move_forward(x, alpha, model_bias)
    x_new = x + 0.5*(x + model_bias + alpha .* x .* abs(x));
end
