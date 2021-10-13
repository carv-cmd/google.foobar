#/bin/bash

show_status () {
		printf '\n\E[38;5;46m>>> Date %s >>> sys:%s >>>\E[00m\n' \
				"$(date +'%Y-%M-%D %T')" \
				"$(vcgencmd measure_temp)"
}

looper () {
		_RANGES=( 10 100 1000 )
		for ix in {2..4}; do
				printf "\n\E[38;5;13m>>> Method Subset: [ ${ix} ]\E[00m\n"
				for inx in "${_RANGES[@]}"; do
						show_status
						local _FOO="$(( "${inx}" * "${1}" ))"
						time python2.7 power-hungry.py -m "${ix}" -i "${_FOO}"
				done
		done
		unset '_RANGES'
}

singles () {
		_RANGES=( 10 100 1000 )
		for inx in "${_RANGES[@]}"; do
				show_status
				local _FOO="$(( "${inx}" * "${2}" ))"
				time python2.7 power-hungry.py -m "${1}" -i "${_FOO}"
		done
		unset '_RANGES'
}


if [ -z "${2}" ]; then
		looper "${1:-1}"
else
		singles "${1:-4}" "${2:-1}"
fi


