#ifndef __VTK_ITK_ITERATIVE_CLOEST_POINT_H__
#define __VTK_ITK_ITERATIVE_CLOEST_POINT_H__

#include <vtkLinearTransform.h>

class vtkITKIterativeCloestPoint : public vtkLinearTransform
{
public:
	static vtkITKIterativeCloestPoint* New();
	vtkTypeMacro(vtkITKIterativeCloestPoint, vtkLinearTransform);
	virtual void PrintSelf(ostream& os, vtkIndent indent) override;


	void SetSource(vtkPoints *points);
	vtkGetObjectMacro(Source, vtkPoints);

	void SetTarget(vtkPoints *points);
	vtkGetObjectMacro(Target, vtkPoints);

	vtkSetMacro(NumberOfIterations, unsigned int);
	vtkGetMacro(NumberOfIterations, unsigned int);

	vtkSetMacro(Mode, int);
	vtkGetMacro(Mode, int);
	void SetModeToRigid() { this->SetMode(RIGID); }
	void SetModeToAffine() { this->SetMode(AFFINE); }

	vtkBooleanMacro(InitializationWithPCA, bool);
	vtkSetMacro(InitializationWithPCA, bool);
	vtkGetMacro(InitializationWithPCA, bool);

	vtkGetMacro(RMS, double);

	enum TRANSFORM_MODE
	{
		RIGID = 0,
		AFFINE,
	};

	virtual void Inverse() override;
	virtual vtkAbstractTransform* MakeTransform();

protected:
	vtkITKIterativeCloestPoint();
	virtual ~vtkITKIterativeCloestPoint() override;

	bool InitializationWithPCA;
	unsigned int NumberOfIterations;
	unsigned int Mode;
	double RMS;

	vtkPoints* Source;
	vtkPoints* Target;

	virtual void InternalDeepCopy(vtkAbstractTransform* transform) override;

	virtual void InternalUpdate() override;

	void ITK_Calculation();
	
private:
	vtkITKIterativeCloestPoint(const vtkITKIterativeCloestPoint&) = delete;
	void operator=(const vtkITKIterativeCloestPoint&) = delete;
	
};

#endif // !__VTK_ITK_ITERATIVE_CLOEST_POINT_H__
