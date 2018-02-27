// itk
#include <itkMesh.h>
#include <itkMeshFileReader.h>
#include <itkMeshFileWriter.h>
#include <itkTransformMeshFilter.h>
#include <itkEuler3DTransform.h>
#include <itkEuclideanDistancePointMetric.h>
#include <itkIdentityTransform.h>
#include <itkEuclideanDistancePointSetToPointSetMetricv4.h>
#include <itkPointSetToPointSetRegistrationMethod.h>
#include <itkLevenbergMarquardtOptimizer.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

#include "vtkITKIterativeCloestPoint.h"
#include "vtkLandmarkTransform.h"

#include <vnl_vector.h>
#include <vnl_matrix.h>
#include <vnl_rotation_matrix.h>
#include <vnl_inverse.h>
#include <vnl_symmetric_eigensystem.h>
#include <vnl_cross.h>

#include <vector>
#include <string>
#include <ctime>
using namespace std;

const unsigned int DIMENSION = 3;
typedef double CoordinateType;
typedef itk::Mesh<CoordinateType, DIMENSION> MeshType;
typedef itk::MeshFileReader<MeshType> MeshFileReader;
typedef itk::MeshFileWriter<MeshType> MeshFileWriter;
typedef itk::EuclideanDistancePointMetric<MeshType, MeshType> EuclideanDistancePointMetric;
typedef itk::EuclideanDistancePointSetToPointSetMetricv4<MeshType, MeshType> EuclideanDistancePointMetricv4;
template <typename TTransform>
using TransformMeshFilter = itk::TransformMeshFilter<MeshType, MeshType, TTransform>;
typedef itk::PointSetToPointSetRegistrationMethod<MeshType, MeshType> PointSetToPointSetRegistration;
typedef itk::Euler3DTransform<PointSetToPointSetRegistration::TransformType::ParametersValueType> Euler3DTransform;
typedef itk::IdentityTransform<EuclideanDistancePointMetric::TransformType::ParametersValueType, DIMENSION> IdentityTransform;
typedef itk::LevenbergMarquardtOptimizer LevenbergMarquardtOptimizer;
//typedef itk::AffineTransform<CoordinateType, DIMENSION> AffineTransform;
/**
 * mesh[0] is the source/moving mesh
 * mesh[1] is the target/fixed mesh
 */
int main(int argc, char** argv)
{
	argv[0] = "../";
	argv[1] = "asdfasdf";
	argv[2] = "asdfasdf";
	string folder("../");
	vector<string> fileNames{
		"Reference with develop tool.vtk",
		"Dental(shift mid point to origin).vtk" };
	//vector<string> fileNames{
	//	"phantom.vtk",
	//	"skin.vtk"
	//};

	vector<vtkSmartPointer<vtkPolyDataReader>> readers;
	for (int i = 0; i < 2; ++i) {
		readers.push_back(vtkSmartPointer<vtkPolyDataReader>::New());
		readers[i]->SetFileName((folder + fileNames[i]).c_str());
		readers[i]->Update();

	}
	vtkSmartPointer<vtkITKIterativeCloestPoint> itk_icp =
		vtkSmartPointer<vtkITKIterativeCloestPoint>::New();
	itk_icp->SetSource(readers[0]->GetOutput()->GetPoints());
	itk_icp->SetTarget(readers[1]->GetOutput()->GetPoints());
	//itk_icp->GetLandmarkTransform()->SetModeToRigidBody();
	//itk_icp->SetMaximumNumberOfIterations(20);
	//itk_icp->StartByMatchingCentroidsOn();
	itk_icp->SetNumberOfIterations(150);
	//itk_icp->InitializationWithPCAOff();
	itk_icp->Update();
	cout << "itk cal RMS: " << itk_icp->GetRMS() << '\n';
	
	cerr << "Finish!\n";
	cin.get();

	return 0;
}

int prevous_main(int argc, char** argv)
{
	argv[0] = "../";
	argv[1] = "asdfasdf";
	argv[2] = "asdfasdf";
	string folder("../");
	vector<string> fileNames{
		"Reference with develop tool.vtk",
		"Dental(shift mid point to origin).vtk" };

	vector<MeshFileReader::Pointer> meshReaders;
	for (int i = 0; i < 2; ++i) {
		meshReaders.push_back(MeshFileReader::New());
	}
	vector<vnl_matrix<CoordinateType>> mesh_matrix(2, vnl_matrix<CoordinateType>());
	vector<vnl_vector_fixed<CoordinateType, DIMENSION>> centers(2, vnl_vector_fixed<CoordinateType, DIMENSION>());
	vector<vnl_vector_fixed<CoordinateType, DIMENSION>> eigenvalues(2, vnl_vector_fixed<CoordinateType, DIMENSION>());
	vector<vnl_matrix_fixed<CoordinateType, DIMENSION, DIMENSION>> eigenvectors(2, vnl_matrix_fixed<CoordinateType, DIMENSION, DIMENSION>());

	for (int i = 0; i < 2; ++i) {
		meshReaders[i]->SetFileName(folder + fileNames[i]);
		meshReaders[i]->Update();
		MeshType::Pointer mesh = meshReaders[i]->GetOutput();
		MeshType::PointsContainer::Pointer points = mesh->GetPoints();
		//vnl_matrix<CoordinateType> matrix(DIMENSION, points->size());
		mesh_matrix[i].set_size(DIMENSION, points->Size());

		for (MeshType::PointsContainer::ConstIterator cit = points->Begin();
			cit != points->End(); ++cit) {
			//mesh_matrix[i].set_column(cit->Index(), cit->Value().GetVnlVector());
			mesh_matrix[i][0][cit->Index()] = cit->Value().GetVnlVector()[0];
			mesh_matrix[i][1][cit->Index()] = cit->Value().GetVnlVector()[1];
			mesh_matrix[i][2][cit->Index()] = cit->Value().GetVnlVector()[2];
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
		covariance /= points->size();



		vnl_symmetric_eigensystem<CoordinateType> eigenSystem(covariance);
		eigenvalues[i] = eigenSystem.D.get_diagonal();
		eigenvectors[i] = eigenSystem.V;
		

	}

	vnl_matrix_fixed< CoordinateType, DIMENSION + 1, DIMENSION + 1> initialization_transform_matrix;
	double initialization_rms = itk::NumericTraits<double>::max();
	MeshType::Pointer initializationMesh;

	for (int i = 0; i < 4; ++i) {
		vector<CoordinateType> opposite(2, 1.f);
		opposite[0] = (i / 2) ? 1.f : -1.f;
		opposite[1] = (i % 2) ? 1.f : -1.f;
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
		vnl_matrix<CoordinateType> source_to_target_matrix =  (transform_matrix * _mesh_matrix[0]).extract(_mesh_matrix[0].rows() - 1, _mesh_matrix[0].columns());
		//vnl_matrix<CoordinateType> source_to_target_matrix =  (_mesh_matrix[0]).extract(_mesh_matrix[0].rows() - 1, _mesh_matrix[0].columns());
		//cout << "source_to_target_matrix" << source_to_target_matrix << '\n';

		MeshType::PointsContainerPointer points = MeshType::PointsContainer::New();
		points->resize(source_to_target_matrix.columns());
		for (MeshType::PointsContainerIterator cit = points->Begin();
			cit != points->End(); ++cit) {
			//cit->Value().GetVnlVector().update(source_to_target_matrix.get_column(cit.Index()));
			cit->Value().GetVnlVector()[0] = source_to_target_matrix[0][cit->Index()];
			cit->Value().GetVnlVector()[1] = source_to_target_matrix[1][cit->Index()];
			cit->Value().GetVnlVector()[2] = source_to_target_matrix[2][cit->Index()];
		}
		MeshType::Pointer sourceToTargetMesh = MeshType::New();
		sourceToTargetMesh->SetPoints(points);
		sourceToTargetMesh->Update();
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
		euclideanDistancev4->SetMovingPointSet(sourceToTargetMesh);
		euclideanDistancev4->SetFixedPointSet(meshReaders[1]->GetOutput());
		double rms = euclideanDistancev4->GetValue();
		if (rms < initialization_rms) {
			initialization_rms = rms;
			initialization_transform_matrix = transform_matrix;
			initializationMesh = sourceToTargetMesh;
		}
		//cout << value << '\n';
		
	//MeshFileWriter::Pointer meshFileWriter =
	//	MeshFileWriter::New();
	////meshFileWriter->SetInput(transformMeshFilter->GetOutput());
	//meshFileWriter->SetInput(sourceToTargetMesh);
	//meshFileWriter->SetFileName(folder + to_string(i) +  "transform.vtk");
	//meshFileWriter->SetFileTypeAsBINARY();
	//meshFileWriter->Write();

	}

	//initialization_transform_matrix.ortho
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
	levenbergMarquardtOptimizer->SetNumberOfIterations(100);
	levenbergMarquardtOptimizer->SetValueTolerance(1e-5);
	levenbergMarquardtOptimizer->SetGradientTolerance(1e-5);
	levenbergMarquardtOptimizer->SetEpsilonFunction(1e-6);

	PointSetToPointSetRegistration::Pointer pointSetToPointSetRegistration =
		PointSetToPointSetRegistration::New();
	pointSetToPointSetRegistration->SetMetric(EuclideanDistancePointMetric::New());
	pointSetToPointSetRegistration->SetOptimizer(levenbergMarquardtOptimizer);
	pointSetToPointSetRegistration->SetMovingPointSet(meshReaders[0]->GetOutput());
	pointSetToPointSetRegistration->SetFixedPointSet(meshReaders[1]->GetOutput());
	pointSetToPointSetRegistration->SetTransform(transform);
	pointSetToPointSetRegistration->SetInitialTransformParameters(transform->GetParameters());
	pointSetToPointSetRegistration->Update();
	cout << "itk_ICP final rms" << levenbergMarquardtOptimizer->GetValue().mean() << '\n';
	//pointSetToPointSetRegistration->

	return 0;
}
