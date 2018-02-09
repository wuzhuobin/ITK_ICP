#ifndef __VTK_PCA_ICP_TRANSFORM_H__
#define __VTK_PCA_ICP_TRANSFORM_H__

#include "vtkIterativeClosestPointTransform.h"

class vtkMatrix4x4;
class vtkPolyData;

class vtkPCA_ICP_Transform : public vtkIterativeClosestPointTransform
{
public:
	static vtkPCA_ICP_Transform* New();
	vtkTypeMacro(vtkPCA_ICP_Transform, vtkIterativeClosestPointTransform);
	virtual void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

	vtkGetMacro(RMS, double);
	vtkGetObjectMacro(TransformPolyData, vtkPolyData);

	vtkGetMacro(ITK_CAL, bool);
	vtkSetMacro(ITK_CAL, bool);
	vtkBooleanMacro(ITK_CAL, bool);

protected:
	vtkPCA_ICP_Transform();
	virtual ~vtkPCA_ICP_Transform() VTK_OVERRIDE;
	virtual void InternalUpdate() VTK_OVERRIDE;
	void ITKCalculate();

	bool ITK_CAL;
	double RMS;
	vtkPolyData* TransformPolyData;

};

#endif // !__VTK_PCA_ICP_TRANSFORM_H__
