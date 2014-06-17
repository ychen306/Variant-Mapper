var variantMapperApp = angular.module('VariantMapper', []);


variantMapperApp.controller('variantMapperController', function($scope){
	$scope.repositories = datasets;
	$scope.readsets = [];
	$scope.pageTokens = [undefined];
	var initRepo = function() {
		for (var repo in datasets) {
			$scope.currentRepo = repo;
			break;
		}
	};

	var initDataset = function() {
		/*
		Select initial Dataset for a newly selected repo.
		*/
		for (var dataset in datasets[$scope.currentRepo]) {
			$scope.currentDataset = dataset;
			break;
		}
	};

	var register = function() {
		var registerUrl = '/register?';
		registerUrl +=  ('repo='+$scope.currentRepo);
		registerUrl += ('&dataset='+$scope.currentDataset);
		registerUrl += ('&readset='+JSON.stringify($scope.currentReadset));
		$.get(registerUrl);
	};

	var loadReadsets = function(pager, oldContext) {
		if (!oldContext || !$scope.loading) {
			getReadsets($scope.currentRepo, $scope.currentDataset, $scope, pager);
		}
	};

	$scope.loadPrevReadsets = function() {
		$scope.readsets = [];
		loadReadsets({last: true});
	};

	$scope.loadNextReadsets = function() {
		$scope.readsets = [];
		loadReadsets({next: true});
	};

	$scope.hasLastPage = function() {
		return $scope.pageTokens.length > 2;
	};

	$scope.hasNextPage = function() {
		return $scope.pageTokens.slice(-1)[0] !== undefined;
	};

	$scope.switchReadset = function(readset) {
		var trimmedReadset = {};
		trimmedReadset.name = readset.name;
		trimmedReadset.id = readset.id;
		$scope.currentReadset = trimmedReadset;
		register();
	};
	$scope.switchDataset = function(dataset) {
		$scope.currentDataset = dataset;
		$scope.readsets = [];
		$scope.pageTokens = [undefined];
		loadReadsets();
	}
	$scope.switchRepo = function(repo) {
		$scope.currentRepo = repo;
		$scope.readsets = [];
		$scope.pageTokens = [undefined];
		initDataset();
		loadReadsets();
	}
	$scope.mapTerm = function(term) {
		$scope.stage = "Downloading reports from Clinvar, NCBI"
		$scope.loadingReports = true;
		mapper = variantMapper($scope.currentRepo, $scope.currentReadset.id, $scope);
		mapper.map(term);
	}
	$scope.noReportsMatched = function() {
		return ($scope.loadingReports === false && (!$scope.matchedReports || $scope.matchedReports.length === 0));
	}
	$scope.canShow = function(report) {
		if ($scope.sigFilter === undefined) 
			return true;
		return ($scope.sigFilter === report.clinicalSignificance)
	}
	$scope.showSigFilter = function() {
		if ($scope.sigFilter === undefined) 
			return "all"
		else
			return $scope.sigFilter;
	}
	$scope.setSigFilter = function(sig) {
		$scope.sigFilter = sig;
	}

	if (!$scope.currentRepo) {
		initRepo();	
	}

	if (!$scope.currentDataset) {
		initDataset();
	}
	loadReadsets();
	$scope.currentReadset = $scope.currentReadset ? $scope.curretnReadset : {name: 'No readset selected'};

});

