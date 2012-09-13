package xmipp.particlepicker.training.gui;

import java.awt.Component;

import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

import xmipp.particlepicker.ParticleCanvas;



public class ParticleCellRenderer implements TableCellRenderer {

	@Override
	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
	{
		ParticleCanvas c = (ParticleCanvas)value;
		return c;
	}
	

	

}
