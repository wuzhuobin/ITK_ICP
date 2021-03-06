// me 
#include "vtkPCA_ICP_Transform.h"

// vtk
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkTable.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPCAStatistics.h>
#include <vtkCenterOfMass.h>
#include <vtkPolyData.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkDoubleArray.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(vtkPCA_ICP_Transform);

vtkPCA_ICP_Transform::vtkPCA_ICP_Transform()
{
	this->TransformPolyData = vtkPolyData::New();
	this->RMS = VTK_DOUBLE_MAX;
	this->UsePCA_Initialization = true;

	this->SourceCenter[0] = 0;
	this->SourceCenter[1] = 0;
	this->SourceCenter[2] = 0;

	this->SourceDirection0[0] = 0;
	this->SourceDirection0[1] = 0;
	this->SourceDirection0[2] = 0;
	this->SourceDirection1[0] = 0;
	this->SourceDirection1[1] = 0;
	this->SourceDirection1[2] = 0;
	this->SourceDirection2[0] = 0;
	this->SourceDirection2[1] = 0;
	this->SourceDirection2[2] = 0;

	this->TargetCenter[0] = 0;
	this->TargetCenter[1] = 0;
	this->TargetCenter[2] = 0;

	this->TargetDirection0[0] = 0;
	this->TargetDirection0[1] = 0;
	this->TargetDirection0[2] = 0;
	this->TargetDirection1[0] = 0;
	this->TargetDirection1[1] = 0;
	this->TargetDirection1[2] = 0;
	this->TargetDirection2[0] = 0;
	this->TargetDirection2[1] = 0;
	this->TargetDirection2[2] = 0;
}

void vtkPCA_ICP_Transform::PrintSelf(ostream & os, vtkIndent indent)
{
}

vtkPCA_ICP_Transform::~vtkPCA_ICP_Transform()
{
	this->TransformPolyData->Delete();
}

void vtkPCA_ICP_Transform::InternalUpdate()
{
	//MarkSegmentation();
	if (this->Source == NULL || !this->Source->GetNumberOfPoints())
	{
		vtkErrorMacro(<< "Can't execute with NULL or empty input");
		return;
	}

	if (!this->Source->IsA("vtkPolyData")) {
		vtkErrorMacro(<< "The Source is not vtkPolyData.");
		return;
	}

	if (this->Target == NULL || !this->Target->GetNumberOfPoints())
	{
		vtkErrorMacro(<< "Can't execute with NULL or empty target");
		return;
	}

	if (!this->Target->IsA("vtkPolyData")) {
		vtkErrorMacro(<< "The Target is not vtkPolyData.");
		return;
	}

	vtkSmartPointer<vtkPolyData> finalPolyData =
		vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkTransform> finalTransform =
		vtkSmartPointer<vtkTransform>::New();
	if (this->UsePCA_Initialization) {
		vtkSmartPointer<vtkCenterOfMass> sourceCenterOfMass = vtkSmartPointer<vtkCenterOfMass>::New();
		sourceCenterOfMass->SetInputData(this->Source);
		sourceCenterOfMass->Update();
		sourceCenterOfMass->GetCenter(this->SourceCenter);

		vtkSmartPointer<vtkCenterOfMass> targetCenterOfMass = vtkSmartPointer<vtkCenterOfMass>::New();
		targetCenterOfMass->SetInputData(this->Target);
		targetCenterOfMass->Update();
		targetCenterOfMass->GetCenter(this->TargetCenter);


		double translation[3] = {
			translation[0] = this->TargetCenter[0] - this->SourceCenter[0],
			translation[1] = this->TargetCenter[1] - this->SourceCenter[1],
			translation[2] = this->TargetCenter[2] - this->SourceCenter[2]
		};
		//cout << "sourceCenter" << " " << sourceCenter[0] << " " << sourceCenter[1] << " " << sourceCenter[2] << " " << endl;

		//double  targetCenter[3];
		//cout << "targetCenter" << " " << targetCenter[0] << " " << targetCenter[1] << " " << targetCenter[2] << " " << endl;

		//double *translation = new double[3];

		//cout << translation[0] << endl;
		//cout << translation[1] << endl;
		//cout << translation[2] << endl;

		vtkSmartPointer<vtkDoubleArray> sourceX = vtkSmartPointer<vtkDoubleArray>::New();
		sourceX->SetNumberOfComponents(1);
		sourceX->SetName("X");

		vtkSmartPointer<vtkDoubleArray> sourceY = vtkSmartPointer<vtkDoubleArray>::New();
		sourceY->SetNumberOfComponents(1);
		sourceY->SetName("Y");

		vtkSmartPointer<vtkDoubleArray> sourceZ = vtkSmartPointer<vtkDoubleArray>::New();
		sourceZ->SetNumberOfComponents(1);
		sourceZ->SetName("Z");

		//cout << "*sourceZ" << " " << *sourceZ << endl;

		//cout << "NumberOfPoints: " << source->GetPoints()->GetNumberOfPoints() << endl;
		for (int i = 0; i < this->Source->GetNumberOfPoints(); ++i) {
			sourceX->InsertNextValue(this->Source->GetPoint(i)[0]);
			sourceY->InsertNextValue(this->Source->GetPoint(i)[1]);
			sourceZ->InsertNextValue(this->Source->GetPoint(i)[2]);
		}
		vtkSmartPointer<vtkTable> sourceTable = vtkSmartPointer<vtkTable>::New();
		sourceTable->AddColumn(sourceX);
		sourceTable->AddColumn(sourceY);
		sourceTable->AddColumn(sourceZ);

		vtkSmartPointer<vtkDoubleArray> sourceEigen[3] = {
			vtkSmartPointer<vtkDoubleArray>::New(),
			vtkSmartPointer<vtkDoubleArray>::New(),
			vtkSmartPointer<vtkDoubleArray>::New() };

		vtkSmartPointer<vtkPCAStatistics> sourcePCA = vtkSmartPointer<vtkPCAStatistics>::New();
		sourcePCA->SetColumnStatus("X", 1);
		sourcePCA->SetColumnStatus("Y", 1);
		sourcePCA->SetColumnStatus("Z", 1);
		sourcePCA->SetInputData(sourceTable);
		sourcePCA->RequestSelectedColumns();
		sourcePCA->SetDeriveOption(true);
		sourcePCA->Update();
		sourcePCA->GetEigenvector(0, sourceEigen[0]);
		sourcePCA->GetEigenvector(1, sourceEigen[1]);
		sourcePCA->GetEigenvector(2, sourceEigen[2]);

		vtkSmartPointer<vtkDoubleArray> targetX = vtkSmartPointer<vtkDoubleArray>::New();
		targetX->SetNumberOfComponents(1);
		targetX->SetName("X");

		vtkSmartPointer<vtkDoubleArray> targetY = vtkSmartPointer<vtkDoubleArray>::New();
		targetY->SetNumberOfComponents(1);
		targetY->SetName("Y");

		vtkSmartPointer<vtkDoubleArray> targetZ = vtkSmartPointer<vtkDoubleArray>::New();
		targetZ->SetNumberOfComponents(1);
		targetZ->SetName("Z");

		//cout << "NumberOfPoints: " << xmlPolyDataReader->GetOutput()->GetPoints()->GetNumberOfPoints() << endl;
		for (int i = 0; i < this->Target->GetNumberOfPoints(); ++i) {
			targetX->InsertNextValue(this->Target->GetPoint(i)[0]);
			targetY->InsertNextValue(this->Target->GetPoint(i)[1]);
			targetZ->InsertNextValue(this->Target->GetPoint(i)[2]);
		}


		vtkSmartPointer<vtkTable> targetTable = vtkSmartPointer<vtkTable>::New();
		targetTable->AddColumn(targetX);
		targetTable->AddColumn(targetY);
		targetTable->AddColumn(targetZ);
		//cout << "targetTable" << *targetTable << endl;

		vtkSmartPointer <vtkDoubleArray> targetEigen[3] = {
			vtkSmartPointer<vtkDoubleArray>::New(),
			vtkSmartPointer<vtkDoubleArray>::New(),
			vtkSmartPointer<vtkDoubleArray>::New() };

		vtkSmartPointer<vtkPCAStatistics> targetPCA = vtkSmartPointer<vtkPCAStatistics>::New();
		targetPCA->SetColumnStatus("X", 1);
		targetPCA->SetColumnStatus("Y", 1);
		targetPCA->SetColumnStatus("Z", 1);
		//targetPCA.SetInputData(0, targetTable);
		targetPCA->SetInputData(targetTable);
		targetPCA->RequestSelectedColumns();
		targetPCA->SetDeriveOption(true);
		targetPCA->Update();
		targetPCA->GetEigenvector(0, targetEigen[0]);
		targetPCA->GetEigenvector(1, targetEigen[1]);
		targetPCA->GetEigenvector(2, targetEigen[2]);


		double distance = VTK_DOUBLE_MAX;//DBL_MAX   	//double distance = Double.MAX_VALUE;

		for (int i = 0; i < 4; i++) {                     //question here ! i< 1 or i<2.
			for (int j = 0; j < 3; j++) {
				this->SourceDirection0[j] = sourceEigen[0]->GetValue(j);
				this->SourceDirection1[j] = sourceEigen[1]->GetValue(j);
				this->SourceDirection2[j] = sourceEigen[2]->GetValue(j);
				this->TargetDirection0[j] = (i % 2 ? -1 : 1) * targetEigen[0]->GetValue(j);
				this->TargetDirection1[j] = (i / 2 ? -1 : 1) * targetEigen[1]->GetValue(j);
				this->TargetDirection2[j] = (i / 2 ? -1 : 1) * targetEigen[2]->GetValue(j);
			}
			//cout << endl;
			//cout << i % 2 << endl;
			//cout << i / 2 << endl;
			//cout << endl;
			//cout << "movingVector" << " " << movingVector[0][0] << " " << movingVector[0][1] << " " << movingVector[0][2] << endl;
			//cout << "movingVector" << " " << movingVector[1][0] << " " << movingVector[1][1] << " " << movingVector[1][2] << endl;
			//cout << "movingValue" << " " << movingValue[0] << endl;
			//cout << "movingValue" << " " << movingValue[1] << endl;
			//cout << "targetVector" << " " << targetVector[0][0] << " " << targetVector[0][1] << " " << targetVector[0][2] << endl;
			//cout << "targetVector" << " " << targetVector[1][0] << " " << targetVector[1][1] << " " << targetVector[1][2] << endl;
			//cout << "targetValue" << " " << targetValue[0] << endl;
			//cout << "targetValue" << " " << targetValue[1] << endl;

			//double *_cross1 = new double[3];
			//double *_cross2 = new double[3];
			// rotate the 0 source direction to target direction, along their cross production vector. 
			// to coincide the 0 direction. 
			double   vector0[3];
			vtkMath::Cross(this->SourceDirection0, this->TargetDirection0, vector0);
			vtkMath::Normalize(vector0);

			//double angle0 = vtkMath::AngleBetweenVectors(this->TargetDirection0, this->SourceDirection0) / vtkMath::Pi()
				//* 180;
			double angle0 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(this->SourceDirection0, this->TargetDirection0));


			// rotate the 1 source direction to 1 target direction, while along the 0 target direction. 
			// in this way, the 0 direction still coincides and 1 direction can coincide as well. 
			double sourceDirection1[4] = { 0,0,0,0 };
			std::copy(this->SourceDirection1, this->SourceDirection1 + 3, sourceDirection1);
			double rotatedSourceDirection1[4];

			vtkSmartPointer<vtkTransform> _transform =
				vtkSmartPointer<vtkTransform>::New();
			_transform->Identity();
			_transform->PostMultiply();
			_transform->RotateWXYZ(angle0, vector0);
			_transform->GetMatrix()->MultiplyPoint(sourceDirection1, rotatedSourceDirection1);




			double   vector1[3];
			//vtkMath::Cross(this->SourceDirection1, this->TargetDirection1, vector1);
			//std::copy(this->TargetDirection0, this->TargetDirection0 + 3, vector1);

			vtkMath::Cross(rotatedSourceDirection1, this->TargetDirection1, vector1);
			double angle1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(rotatedSourceDirection1, this->TargetDirection1));
			//System.out.println("Cross1" + Arrays.toString(_cross1));
			//cout << "_Cross1: " << _cross1[0] << " " << _cross1[1] << " " << _cross1[2] << endl;
			//cout << "_Cross2: " << _cross2[0] << " " << _cross2[1] << " " << _cross2[2] << endl;
			//cout << "Angle1: " << angle1 << endl;
			//cout << "Angle2: " << angle2 << endl;

			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			transform->PostMultiply();

			transform->Translate(-this->SourceCenter[0], -this->SourceCenter[1], -this->SourceCenter[2]);
			transform->RotateWXYZ(angle0, vector0);
			transform->RotateWXYZ(angle1, vector1);
			transform->Translate(this->SourceCenter);
			transform->Translate(translation);
			transform->Update();

			vtkSmartPointer<vtkTransformPolyDataFilter> transformPolyDataFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			transformPolyDataFilter->SetInputData(this->Source);
			transformPolyDataFilter->SetTransform(transform);
			transformPolyDataFilter->Update();

			vtkSmartPointer<vtkDistancePolyDataFilter> distancePolyDataFilter = vtkSmartPointer<vtkDistancePolyDataFilter>::New();
			distancePolyDataFilter->SetInputConnection(0, transformPolyDataFilter->GetOutputPort());
			distancePolyDataFilter->SetInputData(1, this->Target);
			distancePolyDataFilter->Update();

			vtkDataArray *array = distancePolyDataFilter->GetOutput()->GetPointData()->GetScalars();

			double sum = 0;
			//cout << "array->GetNumberOfTuples()" << " " << array->GetNumberOfTuples() << endl;
			for (int id = 0; id < array->GetNumberOfTuples(); ++id) {
				sum += array->GetTuple1(id);
			}

			if (sum < distance) {
				distance = sum;
				finalPolyData = transformPolyDataFilter->GetOutput();
				finalTransform = transform;
			}
		}
		// Lazy way
		vtkDataSet* _source = this->Source;
		this->Source = finalPolyData;
		Superclass::InternalUpdate();
		this->Source = _source;
	}
	else {
		Superclass::InternalUpdate();
	}

	// RMS Calculation;
	finalTransform->Concatenate(this->Matrix);
	finalTransform->GetMatrix(this->Matrix);

	vtkSmartPointer<vtkTransformPolyDataFilter> transformPolyDataFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformPolyDataFilter->SetInputData(this->Source);
	transformPolyDataFilter->SetTransform(finalTransform);
	transformPolyDataFilter->Update();

	vtkSmartPointer<vtkDistancePolyDataFilter> distancePolyDataFilter = vtkSmartPointer<vtkDistancePolyDataFilter>::New();
	distancePolyDataFilter->SetInputConnection(0, transformPolyDataFilter->GetOutputPort());
	distancePolyDataFilter->SetInputData(1, this->Target);
	distancePolyDataFilter->Update();

	vtkDataArray *array = distancePolyDataFilter->GetOutput()->GetPointData()->GetScalars();

	double sum = 0;
	int num = array->GetNumberOfTuples();
	//cout << "array->GetNumberOfTuples()" << " " << array->GetNumberOfTuples() << endl;
	for (int id = 0; id < num; ++id) {
		double d = array->GetTuple1(id);
		d = d*d;
		sum += d;
	}
	sum = sum / num;
	this->RMS = sqrt(sum);
	this->TransformPolyData->ShallowCopy(distancePolyDataFilter->GetOutput());

}
