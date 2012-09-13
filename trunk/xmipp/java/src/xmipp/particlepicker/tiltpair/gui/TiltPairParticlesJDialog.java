package xmipp.particlepicker.tiltpair.gui;

import java.awt.Dimension;
import java.util.List;

import xmipp.particlepicker.ParticleCanvas;
import xmipp.particlepicker.ParticlePickerJFrame;
import xmipp.particlepicker.ParticlesJDialog;
import xmipp.particlepicker.tiltpair.model.UntiltedParticle;
import xmipp.particlepicker.training.model.TrainingParticle;
import xmipp.utils.XmippWindowUtil;
import xmipp.utils.XmippMessage;

public class TiltPairParticlesJDialog extends ParticlesJDialog
{

	public TiltPairParticlesJDialog(ParticlePickerJFrame frame)
	{

		super(frame);

	}

	public void loadParticles(boolean resize)
	{
		int side, rows, columns, width = 0, height = 0;
		List<? extends TrainingParticle> particles = frame.getAvailableParticles();
		side = frame.getSide(frame.getFamily().getSize());

		if (side == 0)
			throw new IllegalArgumentException(XmippMessage.getOutOfBoundsMsg("side"));
		
		if (particles.isEmpty())
		{
			particlespn.removeAll();
			sp.setPreferredSize(new Dimension(200, 200));
			pack();
			return;
		}

		if (resize)
		{
			columns = 1;//Math.min(200, particles.size() * side * 2) / (side * 2);
			rows = (int) Math.ceil(particles.size() / (float) columns);
			width = side * columns * 2;
			height = (side * Math.min(10, rows));
			boolean scroll = (height < rows * side);
			width = width + (scroll ? 40 : 20);
			height = height + (scroll ? 0 : 20);
			sp.setPreferredSize(new Dimension(width, height));
		}
		else
		{
			Dimension d = sp.getSize();
			columns = (int) d.getWidth() / (side * 2);
			rows = (int) Math.ceil((particles.size() / (float) columns));
		}
		particlespn.removeAll();
		particles = frame.getAvailableParticles();
		int index = 0;
		ParticleCanvas c;
		TrainingParticle p;
		UntiltedParticle up;
		columns =  columns * 2;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j+= 2, index ++)
			{
				if (index == particles.size())
					break;
				p = particles.get(index);
				c = p.getParticleCanvas(frame);
				up = (UntiltedParticle) p;
				particlespn.add(c, XmippWindowUtil.getConstraints(constraints, j, i, 1));
				
				if (up.getTiltedParticle() != null)
					particlespn.add(up.getTiltedParticle().getParticleCanvas(frame), XmippWindowUtil.getConstraints(constraints, j + 1, i, 1));
			}
		if (resize)
			pack();
		
	}
}
