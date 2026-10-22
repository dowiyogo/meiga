.PHONY: characterize tolerances forward figures-of-merit

characterize:
	python3 -m characterization.scan_D --d-values 40 80 120 --synthetic-events 5000 --omega-samples 20000

tolerances:
	python3 -m characterization.tolerances.validate

forward:
	python3 -m characterization.forward.synthetic_scenarios

figures-of-merit:
	python3 -m characterization.figures_of_merit.run_fom
