ElectroOpticalSimulator
==========================

The `ElectroOpticalSimulator` is meant to be some sort of calculator that makes our life easy when considering overlap integrals. The current implementation is meant to replicate the coupled mode theory as in :cite:`liu_photonic_2009`. For a given electrical field of the form:

.. math::
	\textbf{E}(\textbf{r}) = \sum_\nu A_\nu (z) \mathbf{\hat{\mathcal{E}}}(x,y) e^{i\beta_\nu z}

we now consider a perturbation in the permitivity tensor of, :math:`\Delta \bar{\mathbf{\epsilon}} (V,x,y)`. We can then formulate a set of coupled mode equations of the form:

.. math::
	\pm \frac{d A(z)}{dz} = \sum_\mu i \kappa _{\nu \mu}A_\mu e^{i(\beta_\mu - \beta_\nu)}

where:

.. math::
	\kappa_{\nu \mu} (V) = \omega \int_{-\infty}^{+\infty} \int_{-\infty}^{+\infty} \hat{\mathcal{E}}_\nu^* \cdot \Delta \bar{\mathbf{\epsilon}} (V,x,y) \cdot \hat{\mathcal{E}}_\mu dxdy

With this :math:`\kappa_{\nu \mu}` we can now retrieve not only the self mode coupling, which translates into a phase and/or amplitude modulation (:math:`\Delta \beta_\nu = \kappa_{\nu\nu}`), but we can also study off diagonal contributions of optical modes that are not purely TE and TM. In fact, we can go a step further and just rotate our permitivity tensor so as to study the influence on the performance of the device by simply considering:

.. math::
	\Delta \bar{\mathbf{\epsilon}} (V,x,y) \to R^{-1}(\theta_x, \theta_y, \theta_z)\Delta \bar{\mathbf{\epsilon}} (V,x,y)R(\theta_x, \theta_y, \theta_z)

where :math:`R^{-1}(\theta_x, \theta_y, \theta_z)` is a generic cartesian rotation matrix.

.. note::
	This package follows the theory of :cite:`liu_photonic_2009` which will be useful for the majority of the cases. However, we do wonder about the limitations of the treatment there. **Could you help us?** In particular we are left with the questions:

	1. The theoretical treatment assumes that the mode fields can be normalized by following an orthonormality relation:
		.. math::
			\iint_{-\infty}^{\infty} 
			\left( 
			\hat{\mathcal{E}}_{\nu} \times \hat{\mathcal{H}}_{\mu}^* 
			+ 
			\hat{\mathcal{E}}_{\mu}^* \times \hat{\mathcal{H}}_{\nu} 
			\right) 
			\cdot \hat{z}\, dx\, dy 
			= \pm \delta_{\nu\mu},

		Will this hold for waveguides with loss as is the general case for doped semiconductor structures?

	2. The coupled mode theory above does not account for lossy modes. Would that change the :math:`\kappa_{\nu \mu}` calculation?
	3. In a lossless waveguide and in the absence of electro-absorption we will have :math:`\Delta \epsilon_{ij}=\Delta \epsilon_{ji}^*`. This is hardly the case in semiconductor structures. Will this jeopardize our calculations? How can we fix this?

.. autoclass:: imodulator.ElectroOpticalSimulator.ElectroOpticalSimulator
   :members:
   :special-members: __init__

References
----------

.. bibliography::
   :filter: docname in docnames
   :style: plain