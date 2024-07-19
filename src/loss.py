import torch
import torch.nn.functional as F

import random
import numpy as np


class BaseLoss:
    eps = 1e-9
    def __init__(self, model):
        self.n_output = len(list(model.clusters[0].parameters())[0])
        self.weight = 1
    @staticmethod
    def compute_distance(is_binary_input, output, target):
        """\
        Compute the distance between target and output with BCE if binary data or MSE for all others.
        """
        if is_binary_input:
            return F.binary_cross_entropy(output, target)
        else:
            return F.mse_loss(output, target)


class SelfEntropyLoss(BaseLoss):
    """
    Entropy regularization to prevent trivial solution.
    """
    def __init__(self, loss_weight=1.0):
        # super().__init__()
        self.prob_layer = torch.nn.Softmax(dim=1)
        self.weight = loss_weight

    def __call__(self, cluster_outputs):
        loss = 0.
        eps = 1e-9
        cluster_outputs = self.prob_layer(cluster_outputs)
        prob_mean = cluster_outputs.mean(dim=0)
        prob_mean[(prob_mean < eps).data] = eps
        # print(prob_mean)
        # print(torch.log(prob_mean))
        loss = -(prob_mean * torch.log(prob_mean)).sum()

        loss *= self.weight
        return loss
    
class DDCLoss(BaseLoss):
    def __init__(self, n_output, loss_weight=1.0):
        # super().__init__()
        self.weight = loss_weight
        self.n_output = n_output
        self.eye = torch.eye(self.n_output, device=torch.device("cuda"))
        self.prob_layer = torch.nn.Softmax(dim=1)

    @staticmethod
    def triu(X):
        """\ 
        Sum of strictly upper triangular part.
        """
        return torch.sum(torch.triu(X, diagonal=1))

    @staticmethod
    def _atleast_epsilon(X, eps=1e-9):
        """
        Ensure that all elements are >= `eps`.
        """
        return torch.where(X < eps, X.new_tensor(eps), X)

    @staticmethod
    def d_cs(A, K, n_clusters):
        """
        Cauchy-Schwarz divergence.
        """
        nom = torch.t(A) @ K @ A
        dnom_squared = torch.unsqueeze(torch.diagonal(nom), -1) @ torch.unsqueeze(
            torch.diagonal(nom), 0
        )

        nom = DDCLoss._atleast_epsilon(nom)
        dnom_squared = DDCLoss._atleast_epsilon(dnom_squared, eps=BaseLoss.eps ** 2)

        d = (
            2
            / (n_clusters * (n_clusters - 1))
            * DDCLoss.triu(nom / torch.sqrt(dnom_squared))
        )
        return d

    @staticmethod
    def kernel_from_distance_matrix(dist, rel_sigma, min_sigma=BaseLoss.eps):
        """\
        Compute a Gaussian kernel matrix from a distance matrix.
        """
        # `dist` can sometimes contain negative values due to floating point errors, so just set these to zero.
        dist = F.relu(dist)
        sigma2 = rel_sigma * torch.median(dist)
        # Disable gradient for sigma
        sigma2 = sigma2.detach()
        sigma2 = torch.where(sigma2 < min_sigma, sigma2.new_tensor(min_sigma), sigma2)
        k = torch.exp(-dist / (2 * sigma2))
        return k

    @staticmethod
    def kernel_from_distance_matrix(dist, rel_sigma, min_sigma=BaseLoss.eps):
        """\
        Compute a Gaussian kernel matrix from a distance matrix.
        """
        # `dist` can sometimes contain negative values due to floating point errors, so just set these to zero.
        dist = F.relu(dist)
        sigma2 = rel_sigma * torch.median(dist)
        # Disable gradient for sigma
        sigma2 = sigma2.detach()
        sigma2 = torch.where(sigma2 < min_sigma, sigma2.new_tensor(min_sigma), sigma2)
        k = torch.exp(-dist / (2 * sigma2))
        return k

    @staticmethod
    def cdist(X, Y):
        """\
        Pairwise distance between rows of X and rows of Y.
        """
        xyT = X @ torch.t(Y)
        x2 = torch.sum(X ** 2, dim=1, keepdim=True)
        y2 = torch.sum(Y ** 2, dim=1, keepdim=True)
        d = x2 - 2 * xyT + torch.t(y2)
        return d

    @staticmethod
    def vector_kernel(x, rel_sigma=0.15):
        """\
        Compute a kernel matrix from the rows of a matrix.
        """
        return DDCLoss.kernel_from_distance_matrix(DDCLoss.cdist(x, x), rel_sigma)

    def __call__(self, hidden, cluster_outputs):
        loss = 0.

        cluster_outputs = self.prob_layer(cluster_outputs)
        hidden_kernel = DDCLoss.vector_kernel(hidden)
        # L_1 loss
        loss = DDCLoss.d_cs(cluster_outputs, hidden_kernel, self.n_output)

        # L_3 loss
        m = torch.exp(-DDCLoss.cdist(cluster_outputs, self.eye))
        loss += DDCLoss.d_cs(m, hidden_kernel, self.n_output)
        loss *= self.weight

        return loss

if __name__ == "__main__":
    sce = SelfEntropyLoss(1.0)
    prob = torch.Tensor([[0.1, 0.9, 0.1], [0.1, 0.9, 0.1]])
    loss = sce(prob)
    print(loss)

    prob = torch.Tensor([[0.6, 0.6, 0.5], [0.4, 0.6, 0.7]])
    loss = sce(prob)
    print(loss)