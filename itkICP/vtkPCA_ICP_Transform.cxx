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

// itk_vnl
#include <vnl_vector.h>
#include <vnl_matrix.h>
#include <vnl_rotation_matrix.h>
#include <vnl_inverse.h>
#include <vnl_symmetric_eigensystem.h>
#include <vnl_cross.h>

// itk
#include <itkPointSet.h>
#include <itkEuler3DTransform.h>
#include <itkEuclideanDistancePointMetric.h>
#include <itkIdentityTransform.h>
#include <itkEuclideanDistancePointSetToPointSetMetricv4.h>
#include <itkPointSetToPointSetRegistrationMethod.h>
#include <itkLevenbergMarquardtOptimizer.h>

// std
#include <vector>


vtkStandardNewMacro(vtkPCA_ICP_Transform);

vtkPCA_ICP_Transform::vtkPCA_ICP_Transform()
{
	this->TransformPolyData = vtkPolyData::New();
	this->ITK_CAL = false;
	this->RMS = VTK_DOUBLE_MAX;
}

void vtkPCA_ICP_Transform::ITKCalculate()
{
	using namespace std;
	typedef double CoordinateType;
	const unsigned int DIMENSION = 3;
	typedef itk::PointSet<CoordinateType, DIMENSION> PointSet;
	typedef itk::EuclideanDistancePointSetToPointSetMetricv4<PointSet, PointSet> EuclideanDistancePointMetricv4;
	typedef itk::EuclideanDistancePointMetric<PointSet, PointSet> EuclideanDistancePointMetric;
	//template <typename TTransform>
	//using TransformMeshFilter = itk::TransformMeshFilter<PointSet, PointSet, TTransform>;
	typedef itk::PointSetToPointSetRegistrationMethod<PointSet, PointSet> PointSetToPointSetRegistration;
	typedef itk::Euler3DTransform<PointSetToPointSetRegistration::TransformType::ParametersValueType> Euler3DTransform;
	typedef itk::IdentityTransform<EuclideanDistancePointMetric::TransformType::ParametersValueType, DIMENSION> IdentityTransform;
	typedef itk::LevenbergMarquardtOptimizer LevenbergMarquardtOptimizer;
	//typedef itk::PointSet<
	vector<vtkDataSet*> dataSets;
	dataSets.resize(2);
	dataSets[0] = this->Source;
	dataSets[1] = this->Target;
	vector<vnl_matrix<CoordinateType>> mesh_matrix(2, vnl_matrix<CoordinateType>());
	vector<vnl_vector_fixed<CoordinateType, DIMENSION>> centers(2, vnl_vector_fixed<CoordinateType, DIMENSION>());
	vector<vnl_vector_fixed<CoordinateType, DIMENSION>> eigenvalues(2, vnl_vector_fixed<CoordinateType, DIMENSION>());
	vector<vnl_matrix_fixed<CoordinateType, DIMENSION, DIMENSION>> eigenvectors(2, vnl_matrix_fixed<CoordinateType, DIMENSION, DIMENSION>());
	for (int i = 0; i < 2; ++i) {
		mesh_matrix[i].set_size(DIMENSION, dataSets[i]->GetNumberOfPoints());

		for (vtkIdType id = 0; id < dataSets[i]->GetNumberOfPoints(); ++id) {
			double* point = dataSets[i]->GetPoint(id);
			mesh_matrix[i][0][id] = point[0];
			mesh_matrix[i][1][id] = point[1];
			mesh_matrix[i][2][id] = point[2];
		}
		vnl_vector<CoordinateType> pointX = mesh_matrix[i].get_row(0);
		vnl_vector<CoordinateType> pointY = mesh_matrix[i].get_row(1);
		vnl_vector<CoordinateType> pointZ = mesh_matrix[i].get_row(2);

		//vnl_vector_fixed<CoordinateType, 3> center(pointX.mean(), pointY.mean(), pointZ.mean());
		centers[i][0] = pointX.mean();
		centers[i][1] = pointY.mean();
		centers[i][2] = pointZ.mean();

		pointX -= pointX.mean();
		pointY -= pointY.mean();
		pointZ -= pointZ.mean();
		vnl_matrix<CoordinateType> mean_matrix(mesh_matrix[i].rows(), mesh_matrix[i].columns());
		mean_matrix.set_row(0, pointX).set_row(1, pointY).set_row(2, pointZ);
		vnl_matrix_fixed<CoordinateType, DIMENSION, DIMENSION> covariance = mean_matrix * (mean_matrix.transpose());
		covariance /= dataSets[i]->GetNumberOfPoints();

		vnl_symmetric_eigensystem<CoordinateType> eigenSystem(covariance);
		eigenvalues[i] = eigenSystem.D.get_diagonal();
		eigenvectors[i] = eigenSystem.V;
	}

	vnl_matrix_fixed< CoordinateType, DIMENSION + 1, DIMENSION + 1> initialization_transform_matrix;
	double initialization_rms = itk::NumericTraits<double>::max();
	PointSet::Pointer initializationPointSet;
	vector<PointSet::Pointer> pointSets;
	pointSets.resize(2);
	for (int i = 0; i < 2; ++i) {
		pointSets[i] = PointSet::New();
		PointSet::PointsContainerPointer points = PointSet::PointsContainer::New();
		points->resize(mesh_matrix[1].columns());
		for (PointSet::PointsContainerIterator cit = points->Begin();
			cit != points->End(); ++cit) {
			//cit->Value().GetVnlVector().update(source_to_target_matrix.get_column(cit.Index()));
			cit->Value().GetVnlVector()[0] = mesh_matrix[i][0][cit->Index()];
			cit->Value().GetVnlVector()[1] = mesh_matrix[i][1][cit->Index()];
			cit->Value().GetVnlVector()[2] = mesh_matrix[i][2][cit->Index()];
		}
		pointSets[i]->SetPoints(points);
	}

	for (int positive_nagative_determinant = 0; positive_nagative_determinant < 4; ++positive_nagative_determinant) {
		vector<CoordinateType> opposite(2, 1.f);
		opposite[0] = (positive_nagative_determinant / 2) ? 1.f : -1.f;
		opposite[1] = (positive_nagative_determinant % 2) ? 1.f : -1.f;
		vector<vnl_vector_fixed<CoordinateType, DIMENSION>> source_axes(2, vnl_vector_fixed<CoordinateType, DIMENSION>());
		source_axes[0] = eigenvectors[0].get_column(2);
		source_axes[1] = eigenvectors[0].get_column(1);
		vector<vnl_vector_fixed<CoordinateType, DIMENSION>> target_axes(2, vnl_vector_fixed<CoordinateType, DIMENSION>());
		target_axes[0] = eigenvectors[1].get_column(2) * opposite[0];
		target_axes[1] = eigenvectors[1].get_column(1) * opposite[1];

		//CoordinateType rotation_angle0 = angle(eigenvectors[0].get_column(2), eigenvectors[1].get_column(2) * opposite[0]);
		CoordinateType rotation_angle0 = angle(source_axes[0], target_axes[0]);
		//vnl_vector_fixed<CoordinateType, DIMENSION> rotation_axis0 = vnl_cross_3d(eigenvectors[0].get_column(2), eigenvectors[1].get_column(2) * opposite[0]).normalize();
		vnl_vector_fixed<CoordinateType, DIMENSION> rotation_axis0 = vnl_cross_3d(source_axes[0], target_axes[0]).normalize();
		//vnl_vector_fixed<double, DIMENSION> double_rotate_axis0;
		//copy(rotation_axis0.begin(), rotation_axis0.data_block() + DIMENSION, double_rotate_axis0.begin());
		rotation_axis0 *= rotation_angle0;
		//cout << "angle0: " << rotation_angle0 << '\n';
		//vnl_matrix_fixed<double, DIMENSION, DIMENSION> double_rotation_matrix0;
		vnl_matrix_fixed<CoordinateType, DIMENSION, DIMENSION> rotation_matrix0;
		vnl_rotation_matrix(rotation_axis0, rotation_matrix0);
		//copy(double_rotation_matrix0.begin(), double_rotation_matrix0.end(), rotation_matrix0.begin());
		//cout << "rotation_matrix0: " << rotation_matrix0 << '\n';


		//vnl_vector_fixed<CoordinateType, 3> rotated_eigenvector = eigenvectors[0].get_column(1);
		//rotated_eigenvector = rotation_matrix0 * rotated_eigenvector;
		//cout << "rotated_eigenvector: " << rotated_eigenvector << '\n';
		source_axes[1] = rotation_matrix0 * source_axes[1];
		//rotated_eigenvector.as_vector() *= rotation_matrix0;

		//CoordinateType rotation_angle1 = angle(rotated_eigenvector, eigenvectors[1].get_column(1) * opposite[1]);
		CoordinateType rotation_angle1 = angle(source_axes[1], target_axes[1]);
		//vnl_vector_fixed<CoordinateType, DIMENSION> rotation_axis1 = vnl_cross_3d(rotated_eigenvector, eigenvectors[1].get_column(1) * opposite[1]).normalize();
		vnl_vector_fixed<CoordinateType, DIMENSION> rotation_axis1 = vnl_cross_3d(source_axes[1], target_axes[1]).normalize();
		//vnl_vector_fixed<double, DIMENSION> double_rotate_axis1;
		//copy(rotation_axis1.begin(), rotation_axis1.end(), double_rotate_axis1.begin());
		rotation_axis1 *= rotation_angle1;

		//vnl_matrix_fixed<double, DIMENSION, DIMENSION> double_rotation_matrix1;
		vnl_matrix_fixed<CoordinateType, DIMENSION, DIMENSION> rotation_matrix1;
		vnl_rotation_matrix(rotation_axis1, rotation_matrix1);
		//copy(double_rotation_matrix1.begin(), double_rotation_matrix1.end(), rotation_matrix1.begin());

		vnl_matrix_fixed<CoordinateType, DIMENSION, DIMENSION> rotation_matrix = rotation_matrix1 * rotation_matrix0;

		vnl_matrix_fixed<CoordinateType, DIMENSION + 1, DIMENSION + 1> _centralization_matrix;
		_centralization_matrix.set_identity();
		for (int i = 0; i < DIMENSION; ++i) {
			_centralization_matrix[i][DIMENSION] = -centers[0][i];

		}

		vnl_matrix_fixed<CoordinateType, DIMENSION + 1, DIMENSION + 1> _centralization_matrix_inverse = vnl_inverse(_centralization_matrix);
		vnl_matrix_fixed<CoordinateType, DIMENSION + 1, DIMENSION + 1> _rotation_matrix;
		_rotation_matrix.set_identity();
		_rotation_matrix.update(rotation_matrix);
		//for (int i = 0; i < 3; ++i) {
		//	for (int j = 0; j < 3; ++j) {
		//		_rotation_matrix[i][j] = rotation_matrix[i][j];
		//	}
		//}

		vnl_matrix_fixed<CoordinateType, DIMENSION + 1, DIMENSION + 1> _translation_matrix;
		vnl_vector_fixed<CoordinateType, DIMENSION> _translation_vector = centers[1] - centers[0];
		_translation_matrix.set_identity();
		for (int i = 0; i < DIMENSION; ++i) {
			_translation_matrix[i][DIMENSION] = _translation_vector[i];
		}
		vnl_matrix_fixed<CoordinateType, DIMENSION + 1, DIMENSION + 1> transform_matrix =
			_translation_matrix *
			_centralization_matrix_inverse *
			_rotation_matrix *
			_centralization_matrix;

		//cout << "transform_matrix: " << transform_matrix << '\n';

		vector<vnl_matrix<CoordinateType>> _mesh_matrix(2, vnl_matrix<CoordinateType>());
		for (int i = 0; i < 1; ++i) {
			_mesh_matrix[i].set_size(mesh_matrix[i].rows() + 1, mesh_matrix[i].columns());
			_mesh_matrix[i].fill(1.0f);
			_mesh_matrix[i].update(mesh_matrix[i]);
		}
		//cout << "_mesh_matrix0: " << _mesh_matrix[0] << '\n';
		//cout << "_mesh_matrix1: " << _mesh_matrix[1] << '\n';
		//cout << transform_matrix * _mesh_matrix[0] << '\n';
		vnl_matrix<CoordinateType> source_to_target_matrix = (transform_matrix * _mesh_matrix[0]).extract(_mesh_matrix[0].rows() - 1, _mesh_matrix[0].columns());
		//vnl_matrix<CoordinateType> source_to_target_matrix =  (_mesh_matrix[0]).extract(_mesh_matrix[0].rows() - 1, _mesh_matrix[0].columns());
		//cout << "source_to_target_matrix" << source_to_target_matrix << '\n';

		PointSet::Pointer sourceToTargetPointSet = PointSet::New();
		PointSet::PointsContainerPointer points = PointSet::PointsContainer::New();
		points->resize(source_to_target_matrix.columns());
		for (PointSet::PointsContainerIterator cit = points->Begin();
			cit != points->End(); ++cit) {
			//cit->Value().GetVnlVector().update(source_to_target_matrix.get_column(cit.Index()));
			cit->Value().GetVnlVector()[0] = source_to_target_matrix[0][cit->Index()];
			cit->Value().GetVnlVector()[1] = source_to_target_matrix[1][cit->Index()];
			cit->Value().GetVnlVector()[2] = source_to_target_matrix[2][cit->Index()];
		}
		sourceToTargetPointSet->SetPoints(points);
		

		//cout << *sourceToTargetMesh << endl;
		//cout << "rotation_matrix: " << rotation_matrix << '\n';


		//Euler3DTransform::MatrixType rotationMatrix;
		//rotationMatrix.GetVnlMatrix() = rotation_matrix;
		//Euler3DTransform::Pointer euler3DTransform = Euler3DTransform::New();
		//euler3DTransform->SetCenter(centers[0].data_block());
		//euler3DTransform->SetMatrix(rotationMatrix);
		//euler3DTransform->SetTranslation((centers[1] - centers[0]).data_block());
		////cout << *euler3DTransform << '\n';

		//TransformMeshFilter<Euler3DTransform>::Pointer transformMeshFilter =
		//	TransformMeshFilter<Euler3DTransform>::New();
		//transformMeshFilter->SetInput(meshReaders[0]->GetOutput());
		//transformMeshFilter->SetTransform(euler3DTransform);
		//transformMeshFilter->Update();

		EuclideanDistancePointMetricv4::Pointer euclideanDistancev4 = EuclideanDistancePointMetricv4::New();
		euclideanDistancev4->SetMovingPointSet(sourceToTargetPointSet);
		euclideanDistancev4->SetFixedPointSet(pointSets[1]);
		double rms = euclideanDistancev4->GetValue();
		//cout << rms << '\n';
		//cout << "time used: " << 1000.0 * (clock() - start) / CLOCKS_PER_SEC << '\n';
		if (rms < initialization_rms) {
			initialization_rms = rms;
			initialization_transform_matrix = transform_matrix;
			initializationPointSet = sourceToTargetPointSet;
		}
	}


	vnl_matrix_fixed<Euler3DTransform::ParametersValueType, DIMENSION + 1, DIMENSION + 1> double_initialization_transform_matrix;
	copy(initialization_transform_matrix.begin(), initialization_transform_matrix.end(), double_initialization_transform_matrix.begin());

	Euler3DTransform::OffsetType offset;
	offset.SetVnlVector(double_initialization_transform_matrix.get_column(3));
	Euler3DTransform::MatrixType matrix;
	matrix.GetVnlMatrix().update(double_initialization_transform_matrix.extract(DIMENSION, DIMENSION));
	//matrix.GetVnlMatrix()
	Euler3DTransform::Pointer transform = Euler3DTransform::New();
	transform->SetOffset(offset);
	transform->SetMatrix(matrix);

	LevenbergMarquardtOptimizer::ScalesType scales;
	//(transform->GetNumberOfParameters());
	scales.SetSize(transform->GetNumberOfParameters());
	scales.Fill(0.01);


	LevenbergMarquardtOptimizer::Pointer levenbergMarquardtOptimizer =
		LevenbergMarquardtOptimizer::New();
	levenbergMarquardtOptimizer->SetUseCostFunctionGradient(false);
	levenbergMarquardtOptimizer->SetScales(scales);
	levenbergMarquardtOptimizer->SetNumberOfIterations(this->NumberOfIterations);
	levenbergMarquardtOptimizer->SetValueTolerance(1e-5);
	levenbergMarquardtOptimizer->SetGradientTolerance(1e-5);
	levenbergMarquardtOptimizer->SetEpsilonFunction(1e-6);

	PointSetToPointSetRegistration::Pointer pointSetToPointSetRegistration =
		PointSetToPointSetRegistration::New();
	pointSetToPointSetRegistration->SetMetric(EuclideanDistancePointMetric::New());
	pointSetToPointSetRegistration->SetOptimizer(levenbergMarquardtOptimizer);
	pointSetToPointSetRegistration->SetMovingPointSet(pointSets[0]);
	pointSetToPointSetRegistration->SetFixedPointSet(pointSets[1]);
	pointSetToPointSetRegistration->SetTransform(transform);
	pointSetToPointSetRegistration->SetInitialTransformParameters(transform->GetParameters());
	pointSetToPointSetRegistration->Update();
	this->RMS = levenbergMarquardtOptimizer->GetValue().mean();

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 4; ++j) {
			if (j < 3) {
				this->Matrix->SetElement(i, j, transform->GetMatrix()[i][j]);
			}
			else {
				this->Matrix->SetElement(i, j, transform->GetOffset()[i]);
			}
		}
	}
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

	if (this->Target == NULL || !this->Target->GetNumberOfPoints())
	{
		vtkErrorMacro(<< "Can't execute with NULL or empty target");
		return;
	}

	if (this->ITK_CAL) {
		this->ITKCalculate();
		return;
	}


	if (!this->Source->IsA("vtkPolyData")) {
		vtkErrorMacro(<< "The Source is not vtkPolyData.");
		return;
	}

	if (!this->Target->IsA("vtkPolyData")) {
		vtkErrorMacro(<< "The Target is not vtkPolyData.");
		return;
	}

	double SourceCenter[3];
	double TargetCenter[3];
	
	vtkSmartPointer<vtkCenterOfMass> sourceCenterOfMass = vtkSmartPointer<vtkCenterOfMass>::New();
	sourceCenterOfMass->SetInputData(this->Source);
	sourceCenterOfMass->Update();
	sourceCenterOfMass->GetCenter(SourceCenter);

	vtkSmartPointer<vtkCenterOfMass> targetCenterOfMass = vtkSmartPointer<vtkCenterOfMass>::New();
	targetCenterOfMass->SetInputData(this->Target);
	targetCenterOfMass->Update();
	targetCenterOfMass->GetCenter(TargetCenter);


	double translation[3] = {
		translation[0] = TargetCenter[0] - SourceCenter[0],
		translation[1] = TargetCenter[1] - SourceCenter[1],
		translation[2] = TargetCenter[2] - SourceCenter[2]
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
	//sourcePCA.SetInputData(0, movingTable);
	sourcePCA->SetInputData(sourceTable);
	sourcePCA->RequestSelectedColumns();
	sourcePCA->SetDeriveOption(true);
	sourcePCA->Update();
	sourcePCA->GetEigenvector(0, sourceEigen[0]);
	sourcePCA->GetEigenvector(1, sourceEigen[1]);
	sourcePCA->GetEigenvector(2, sourceEigen[2]);
	//sourcePCA->GetEigenvector(1, movingEigen[1]);

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
	//vtkSmartPointer<vtkDoubleArray> targetEigen2 = vtkSmartPointer<vtkDoubleArray>::New();
	//vtkSmartPointer<vtkDoubleArray> targetEigen[2] = { targetEigen1,targetEigen2 };

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
	//targetPCA->GetEigenvector(1, targetEigen[1]);


	//const unsigned int N = 3;//2 ROWS 3 COLOMNS
	//double this->SourceDirection[N];
	//double this->TargetDirection[N];


	//double movingValue[2];
	//double targetValue[2];
	//movingValue[0] = sourcePCA->GetEigenvalue(0);
	////movingValue[1] = sourcePCA->GetEigenvalue(1);
	////targetValue[1] = targetPCA->GetEigenvalue(1);
	//targetValue[0] = targetPCA->GetEigenvalue(0);
 
	vtkSmartPointer<vtkPolyData> finalPolyData;
	//= vtkSmartPointer<vtkPolyData>::New();// define initializationPolyData = null
	vtkSmartPointer<vtkTransform> finalTransform;
	//= vtkSmartPointer<vtkTransform>::New();

	double distance = VTK_DOUBLE_MAX;//DBL_MAX   	//double distance = Double.MAX_VALUE;
	double SourceDirection0[3];
	double SourceDirection1[3];
	double TargetDirection0[3];
	double TargetDirection1[3];
	for (int i = 0; i < 4; i++) {                     //question here ! i< 1 or i<2.
		for (int j = 0; j < 3; j++) {
			SourceDirection0[j] = sourceEigen[0]->GetValue(j);
			SourceDirection1[j] = sourceEigen[1]->GetValue(j);
			TargetDirection0[j] = (i % 2 ? -1 : 1) * targetEigen[0]->GetValue(j);
			TargetDirection1[j] = (i / 2 ? -1 : 1) * targetEigen[1]->GetValue(j);
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
		vtkMath::Cross(SourceDirection0, TargetDirection0, vector0);
		vtkMath::Normalize(vector0);

		//double angle0 = vtkMath::AngleBetweenVectors(this->TargetDirection0, this->SourceDirection0) / vtkMath::Pi()
			//* 180;
		double angle0 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(SourceDirection0, TargetDirection0));


		// rotate the 1 source direction to 1 target direction, while along the 0 target direction. 
		// in this way, the 0 direction still coincides and 1 direction can coincide as well. 
		double sourceDirection1[4] = { 0,0,0,0 };
		std::copy(SourceDirection1, SourceDirection1 + 3, sourceDirection1);
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
		
		vtkMath::Cross(rotatedSourceDirection1, TargetDirection1, vector1);
		double angle1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(rotatedSourceDirection1, TargetDirection1));
		//System.out.println("Cross1" + Arrays.toString(_cross1));
		//cout << "_Cross1: " << _cross1[0] << " " << _cross1[1] << " " << _cross1[2] << endl;
		//cout << "_Cross2: " << _cross2[0] << " " << _cross2[1] << " " << _cross2[2] << endl;
		//cout << "Angle1: " << angle1 << endl;
		//cout << "Angle2: " << angle2 << endl;

		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		transform->PostMultiply();
		
		transform->Translate(-SourceCenter[0], -SourceCenter[1], -SourceCenter[2]);
		transform->RotateWXYZ(angle0, vector0);
		transform->RotateWXYZ(angle1, vector1);
		transform->Translate(SourceCenter);
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

	//vtkSmartPointer<vtkMatrix4x4> matrix =
	//	vtkSmartPointer<vtkMatrix4x4>::New();
	//matrix->DeepCopy(this->Matrix);

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
